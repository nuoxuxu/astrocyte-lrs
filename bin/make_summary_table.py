#!/usr/bin/env python3
import argparse
import csv
import re
import sys

import polars as pl


def parse_gtf_attributes(attr_str):
    attrs = {}
    for match in re.finditer(r'(\w+) "([^"]*)"', attr_str):
        attrs[match.group(1)] = match.group(2)
    return attrs


def get_attr(attr_str, key):
    """Extract a single GTF attribute value by key via an anchored search.

    The collaborator GTF contains malformed empty gene_name attributes written
    with an unbalanced quote (`gene_name ";`). That shifts the general
    parse_gtf_attributes regex and causes the following attributes (e.g. ORF_id)
    to be dropped. Anchoring on the literal key avoids that misalignment.
    """
    match = re.search(re.escape(key) + r' "([^"]*)"', attr_str)
    return match.group(1) if match else None


def parse_orf_locations(gtf_path):
    """Return ORF_id -> (seqname, min_pos, max_pos) from start_codon/stop_codon features."""
    locations = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            if fields[2] not in ('start_codon', 'stop_codon'):
                continue
            seqname = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            orf_id = get_attr(fields[8], 'ORF_id')
            if not orf_id:
                continue
            if orf_id not in locations:
                locations[orf_id] = (seqname, start, end)
            else:
                _, cur_min, cur_max = locations[orf_id]
                locations[orf_id] = (seqname, min(cur_min, start), max(cur_max, end))
    return locations


def parse_phylocsf_scores(gtf_path):
    """Return ORF_id -> (score, power) from transcript-level features.

    The collaborator GTF (input to phylocsfpp) uses ORF_id as transcript_id,
    so we key by transcript_id here which equals ORF_id in the ribotie CSV.
    """
    scores = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'transcript':
                continue
            attrs = parse_gtf_attributes(fields[8])
            orf_id = attrs.get('transcript_id')
            if orf_id:
                scores[orf_id] = (
                    attrs.get('phylocsf_score_weighted_mean', ''),
                    attrs.get('phylocsf_power_mean', ''),
                )
    return scores


def parse_coloc_gene_names(tsv_path, exclude_disease='Height'):
    """Return set of gene symbols from an sQTL-coloc matches TSV.

    Rows whose 'disease' column equals exclude_disease are dropped first.
    The 'gene_names' column is semicolon-separated (and may mix Ensembl IDs
    with symbols), so each token is added individually.
    """
    genes = set()
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if exclude_disease is not None and row.get('disease') == exclude_disease:
                continue
            for token in row.get('gene_names', '').split(';'):
                token = token.strip()
                if token:
                    genes.add(token)
    return genes


def parse_concordant_ribo(concordant_csv_path):
    """Return id -> ribo_log2FC from the collaborator concordant CSV.

    The collaborator file keys ORFs by its second column ('id'), which matches
    the ORF_id (first column) of the summary table. This file is pre-filtered to
    ribo-significant ORFs (ribo_sig == True for every row), so an ORF's presence
    in the returned dict means it has a significant Ribo-seq change, and its
    ribo_log2FC gives the direction of that change.
    """
    values = {}
    with open(concordant_csv_path) as f:
        for row in csv.DictReader(f):
            try:
                values[row['id']] = float(row['ribo_log2FC'])
            except (ValueError, TypeError):
                continue
    return values


def classify_riboseq(ribo_fc, tx_fc, tx_q, q_threshold=0.05):
    """Classify an ORF's transcriptomic dynamics relative to Ribo-seq.

    ribo_fc: ribo_log2FC if the ORF is Ribo-seq significant, else None.
    tx_fc, tx_q: a transcriptomic effect size and its q-value, both from the
    isoform features table; None when unavailable. Called once with the isoform
    expression change (iso_log2_fold_change / iso_q_value) to build
    `riboseq_iso_q`, and once with the isoform usage/switch change
    (dIF / isoform_switch_q_value) to build `riboseq_switch_q`. These columns
    correspond to the iso_q and switch_q columns of the collaborator
    concordant CSV (validated to equal the features-table values).

    - 'concordant': transcript change significant (tx_q < threshold) and Ribo-seq
      significant in the same direction.
    - 'discordant': transcript change significant and Ribo-seq significant in the
      opposite direction.
    - 'buffering':  transcript change significant but no significant Ribo-seq change.
    - 'other':      everything else.
    """
    tx_sig = tx_q is not None and tx_q < q_threshold
    ribo_sig = ribo_fc is not None
    if tx_sig and ribo_sig and tx_fc is not None and tx_fc != 0 and ribo_fc != 0:
        return 'concordant' if (ribo_fc > 0) == (tx_fc > 0) else 'discordant'
    if tx_sig and not ribo_sig:
        return 'buffering'
    return 'other'


def parse_pfam_domain_hits(pfam_csv_path, evalue_threshold=1.0):
    """Return set of ORF_ids (pfam seq_id) with at least one Pfam hit below evalue_threshold.

    ORFs absent from the pfam CSV are simply not in the returned set, so callers
    treat them as having no domain (Domain = 0).
    """
    hits = set()
    with open(pfam_csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if float(row['evalue']) < evalue_threshold:
                hits.add(row['seq_id'])
    return hits


def _cds_signature(chrom, strand, intervals):
    """Canonical signature for a CDS: (chrom, strand, sorted unique (start, end) pairs)."""
    return (chrom, strand, tuple(sorted(set(intervals))))


def _splice_junctions(chrom, exons):
    """Set of (chrom, donor, acceptor) splice junctions for a transcript's exons.

    Each junction is the intron between two consecutive exons (sorted by genomic
    position): donor = upstream exon end, acceptor = downstream exon start. This
    matches the intron convention used by ggtranscript::to_intron in the plotting
    script, so the two stay consistent.
    """
    segs = sorted(exons)
    return {(chrom, segs[i][1], segs[i + 1][0]) for i in range(len(segs) - 1)}


def load_gencode_data(gtf_path):
    """Parse GENCODE GTF in one pass.

    Returns:
        cds_signatures: set of (chrom, strand, sorted unique (start, end) pairs)
        gene_biotypes:  dict  gene_id  -> gene_type
        cds_by_tx:      dict  transcript_id -> {'strand': str, 'segs': [(start, end)]}
        gene_names:     dict  gene_id  -> gene_name
        junctions:      set of (chrom, donor, acceptor) splice junctions across
                        all GENCODE transcripts (every biotype)
    """
    tx_cds = {}
    tx_exons = {}
    gene_biotypes = {}
    gene_names = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            feat = fields[2]
            if feat == 'gene':
                attrs = parse_gtf_attributes(fields[8])
                g = attrs.get('gene_id')
                t = attrs.get('gene_type')
                n = attrs.get('gene_name')
                if g and t:
                    gene_biotypes[g] = t
                if g and n:
                    gene_names[g] = n
            elif feat == 'CDS':
                chrom, strand = fields[0], fields[6]
                start, end = int(fields[3]), int(fields[4])
                attrs = parse_gtf_attributes(fields[8])
                tx_id = attrs.get('transcript_id')
                if not tx_id:
                    continue
                if tx_id not in tx_cds:
                    tx_cds[tx_id] = {'chrom': chrom, 'strand': strand, 'segs': [], 'intervals': []}
                tx_cds[tx_id]['segs'].append((start, end))
                tx_cds[tx_id]['intervals'].append((start, end))
            elif feat == 'exon':
                chrom = fields[0]
                start, end = int(fields[3]), int(fields[4])
                tx_id = get_attr(fields[8], 'transcript_id')
                if not tx_id:
                    continue
                if tx_id not in tx_exons:
                    tx_exons[tx_id] = {'chrom': chrom, 'segs': []}
                tx_exons[tx_id]['segs'].append((start, end))

    cds_signatures = {
        _cds_signature(v['chrom'], v['strand'], v['intervals'])
        for v in tx_cds.values()
    }
    cds_by_tx = {tx: {'strand': v['strand'], 'segs': v['segs']} for tx, v in tx_cds.items()}
    junctions = set()
    for v in tx_exons.values():
        junctions |= _splice_junctions(v['chrom'], v['segs'])
    return cds_signatures, gene_biotypes, cds_by_tx, gene_names, junctions


def parse_study2_cds_signatures(gtf_path):
    """Return set of CDS signatures from study2 ORF GTF (uses lowercase orf_id attribute)."""
    orf_cds = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#') or '\t' not in line:
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            chrom, strand = fields[0], fields[6]
            start, end = int(fields[3]), int(fields[4])
            attrs = parse_gtf_attributes(fields[8])
            orf_id = attrs.get('orf_id')
            if not orf_id:
                continue
            if orf_id not in orf_cds:
                orf_cds[orf_id] = {'chrom': chrom, 'strand': strand, 'intervals': []}
            orf_cds[orf_id]['intervals'].append((start, end))
    return {_cds_signature(v['chrom'], v['strand'], v['intervals']) for v in orf_cds.values()}


def parse_collaborator_cds_signatures(gtf_path):
    """Return ORF_id -> CDS signature from the collaborator GTF."""
    orf_cds = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            chrom, strand = fields[0], fields[6]
            start, end = int(fields[3]), int(fields[4])
            orf_id = get_attr(fields[8], 'ORF_id')
            if not orf_id:
                continue
            if orf_id not in orf_cds:
                orf_cds[orf_id] = {'chrom': chrom, 'strand': strand, 'intervals': []}
            orf_cds[orf_id]['intervals'].append((start, end))
    return {
        orf_id: _cds_signature(v['chrom'], v['strand'], v['intervals'])
        for orf_id, v in orf_cds.items()
    }


def parse_collaborator_junctions(gtf_path):
    """Return ORF_id -> set of (chrom, donor, acceptor) splice junctions.

    Built from the collaborator transcript's exon chain (in this GTF each ORF is
    a transcript whose transcript_id equals its ORF_id). transcript_id is read
    rather than ORF_id because it precedes the occasionally malformed gene_name
    attribute and so is always parseable.
    """
    tx_exons = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            chrom = fields[0]
            start, end = int(fields[3]), int(fields[4])
            orf_id = get_attr(fields[8], 'transcript_id')
            if not orf_id:
                continue
            if orf_id not in tx_exons:
                tx_exons[orf_id] = {'chrom': chrom, 'segs': []}
            tx_exons[orf_id]['segs'].append((start, end))
    return {
        orf_id: _splice_junctions(v['chrom'], v['segs'])
        for orf_id, v in tx_exons.items()
    }


def load_orfanage_data(gtf_path):
    """Parse orfanage GTF in one pass.

    Returns:
        canonical_positions: dict  transcript_id -> (canonical_TIS_idx, canonical_LTS_idx, canonical_TTS_idx)
        iso_exons:           dict  transcript_id -> {'strand': str, 'exons': [(start, end)]}

    Implements the coordinate logic from transcript_transformer/data.py:306-352.
    CDS features in orfanage GTF do not include the stop codon.
    """
    exons   = {}
    cds_segs = {}
    strands  = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            feat = fields[2]
            if feat not in ('exon', 'CDS'):
                continue
            strand = fields[6]
            start, end = int(fields[3]), int(fields[4])
            attrs = parse_gtf_attributes(fields[8])
            tx_id = attrs.get('transcript_id')
            if not tx_id:
                continue
            strands[tx_id] = strand
            if feat == 'exon':
                exons.setdefault(tx_id, []).append((start, end))
            else:
                cds_segs.setdefault(tx_id, []).append((start, end))

    # build iso_exons for GENCODE projection
    iso_exons = {
        tx: {'strand': strands[tx], 'exons': exons[tx]}
        for tx in exons
    }

    # build canonical positions from orfanage CDS
    canonical_positions = {}
    for tx_id, cds_list in cds_segs.items():
        strand    = strands.get(tx_id, '+')
        exon_list = exons.get(tx_id, [])
        if not exon_list:
            continue
        if strand == '+':
            exon_sorted = sorted(exon_list, key=lambda x: x[0])
            tis = min(s for s, e in cds_list)
        else:
            exon_sorted = sorted(exon_list, key=lambda x: -x[1])
            tis = max(e for s, e in cds_list)
        cum, tis_idx = 0, -1
        for es, ee in exon_sorted:
            if strand == '+' and es <= tis <= ee:
                tis_idx = cum + (tis - es); break
            elif strand == '-' and es <= tis <= ee:
                tis_idx = cum + (ee - tis); break
            cum += ee - es + 1
        if tis_idx == -1:
            continue
        total_cds_len = sum(e - s + 1 for s, e in cds_list)
        tts_idx = tis_idx + total_cds_len
        canonical_positions[tx_id] = (tis_idx, tts_idx - 1, tts_idx)

    return canonical_positions, iso_exons


def project_gencode_cds(iso_exon_data, cds_data):
    """Project a GENCODE CDS onto an Iso-Seq transcript's coordinate system.

    Returns (canonical_TIS_idx, canonical_LTS_idx, canonical_TTS_idx) or None
    if the GENCODE TIS does not fall within any Iso-Seq exon.
    """
    strand   = iso_exon_data['strand']
    exons    = iso_exon_data['exons']
    cds_segs = cds_data['segs']

    if strand == '+':
        exon_sorted = sorted(exons, key=lambda x: x[0])
        tis = min(s for s, e in cds_segs)
    else:
        exon_sorted = sorted(exons, key=lambda x: -x[1])
        tis = max(e for s, e in cds_segs)

    cum, tis_idx = 0, -1
    for es, ee in exon_sorted:
        if strand == '+' and es <= tis <= ee:
            tis_idx = cum + (tis - es); break
        elif strand == '-' and es <= tis <= ee:
            tis_idx = cum + (ee - tis); break
        cum += ee - es + 1

    if tis_idx == -1:
        return None

    # CDS length = intersection of GENCODE CDS segments with Iso-Seq exons
    cds_len = 0
    for cs, ce in cds_segs:
        for es, ee in exons:
            ovl_s, ovl_e = max(cs, es), min(ce, ee)
            if ovl_s <= ovl_e:
                cds_len += ovl_e - ovl_s + 1

    tts_idx = tis_idx + cds_len
    return (tis_idx, tts_idx - 1, tts_idx)


def classify_orf_type(TIS_idx, TTS_idx, canonical, tx_id, gene_for_tx, gene_biotypes):
    """Classify ORF type using the transcript_transformer decision tree (processing.py:607-641).

    canonical: (canonical_TIS_idx, canonical_LTS_idx, canonical_TTS_idx) or None.
    Falls back to lncRNA-ORF / varRNA-ORF using gene biotype when canonical is None.
    """
    LTS_idx = TTS_idx - 1
    if canonical is None:
        gene_id = gene_for_tx.get(tx_id)
        biotype = gene_biotypes.get(gene_id, '') if gene_id else ''
        return 'lncRNA-ORF' if biotype == 'lncRNA' else 'varRNA-ORF'
    canonical_TIS_idx, canonical_LTS_idx, canonical_TTS_idx = canonical
    if canonical_TIS_idx == TIS_idx:
        if canonical_LTS_idx == LTS_idx:
            return 'annotated CDS'
        elif canonical_TTS_idx < TTS_idx:
            return 'C-terminal extension'
        else:
            return 'C-terminal truncation'
    elif canonical_TTS_idx < TIS_idx:
        return 'dORF'
    elif canonical_TIS_idx >= TTS_idx:
        return 'uORF'
    elif canonical_TIS_idx > TIS_idx:
        if canonical_TTS_idx == TTS_idx:
            return 'N-terminal extension'
        else:
            return 'uoORF'
    elif canonical_TTS_idx < TTS_idx:
        return 'doORF'
    else:
        if canonical_TTS_idx == TTS_idx:
            return 'N-terminal truncation'
        else:
            return 'intORF'


def main():
    parser = argparse.ArgumentParser(description='Build per-ORF summary table.')
    parser.add_argument('ribotie_csv', help='RiboTIE CSV (ribotie_cpm1_3sample.csv)')
    parser.add_argument('collaborator_gtf', help='Collaborator GTF with ORF_id attributes')
    parser.add_argument('phylocsfpp_gtf', help='PhyloCSF++ annotated GTF')
    parser.add_argument('pfam_csv', help='Merged Pfam scan results CSV')
    parser.add_argument('annotation_gtf', help='GENCODE annotation GTF')
    parser.add_argument('orfanage_gtf', help='Orfanage GTF for canonical CDS reference')
    parser.add_argument('final_classification', help='SQANTI3 final_classification.parquet for lncRNA biotype lookup')
    parser.add_argument('isoform_features_csv', help='IsoformSwitchAnalyzeR isoformFeatures.csv')
    parser.add_argument('study2_gtf', help='Study 2 reference ORF GTF (data/study2_orfs.gtf)')
    parser.add_argument('leafcutter_coloc', help='leafcutter_sqtl_coloc_matches.tsv')
    parser.add_argument('novel_coding_junction_coloc', help='novel_coding_junction_sqtl_coloc_matches.tsv')
    parser.add_argument('concordant_csv', help='Collaborator concordant_all_three_exp.csv (id, ribo_log2FC; ribo-significant ORFs)')
    parser.add_argument('-o', '--output', default='summary_table.tsv')
    args = parser.parse_args()

    orf_locations    = parse_orf_locations(args.collaborator_gtf)
    tx_scores        = parse_phylocsf_scores(args.phylocsfpp_gtf)
    pfam_hits        = parse_pfam_domain_hits(args.pfam_csv)
    study2_signatures = parse_study2_cds_signatures(args.study2_gtf)
    leafcutter_coloc_genes = parse_coloc_gene_names(args.leafcutter_coloc)
    novel_junction_coloc_genes = parse_coloc_gene_names(args.novel_coding_junction_coloc)
    ribo_log2fc = parse_concordant_ribo(args.concordant_csv)

    (gencode_signatures, gene_biotypes, gencode_cds_by_tx, gene_names,
     gencode_junctions) = load_gencode_data(args.annotation_gtf)
    collaborator_signatures = parse_collaborator_cds_signatures(args.collaborator_gtf)
    collaborator_junctions = parse_collaborator_junctions(args.collaborator_gtf)
    canonical_positions, iso_exons = load_orfanage_data(args.orfanage_gtf)

    clf = pl.read_parquet(args.final_classification).select(
        ['isoform', 'associated_gene', 'associated_transcript']
    )
    gene_for_tx          = dict(zip(clf['isoform'], clf['associated_gene']))
    isoform_to_gencode_tx = dict(zip(clf['isoform'], clf['associated_transcript']))

    switch_q = {}
    gene_q = {}
    if_overall = {}
    iso_log2fc = {}
    iso_q_value = {}
    dif = {}
    switch_q_value = {}
    with open(args.isoform_features_csv) as f:
        for row in csv.DictReader(f):
            iso = row['isoform_id']
            try:
                switch_q[iso] = 1 if float(row.get('isoform_switch_q_value', 'NA')) < 0.05 else 0
            except (ValueError, TypeError):
                switch_q[iso] = 0
            try:
                switch_q_value[iso] = float(row.get('isoform_switch_q_value', 'NA'))
            except (ValueError, TypeError):
                pass
            try:
                dif[iso] = float(row.get('dIF', 'NA'))
            except (ValueError, TypeError):
                pass
            try:
                gene_q[iso] = 1 if float(row.get('gene_q_value', 'NA')) > 0.05 else 0
            except (ValueError, TypeError):
                gene_q[iso] = 0
            try:
                if_overall[iso] = 1 if float(row.get('IF_overall', 'NA')) > 0.1 else 0
            except (ValueError, TypeError):
                if_overall[iso] = 0
            try:
                iso_log2fc[iso] = float(row.get('iso_log2_fold_change', 'NA'))
            except (ValueError, TypeError):
                pass
            try:
                iso_q_value[iso] = float(row.get('iso_q_value', 'NA'))
            except (ValueError, TypeError):
                pass

    with open(args.ribotie_csv) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in)
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow([
            'ORF_id', 'gene_name', 'location', 'transcript_id', 'transcript_len',
            'start_codon', 'stop_codon', 'strand',
            'ORF_type_ORFanage', 'ORF_type_RiboTIE',
            'phylocsf_score_weighted_mean', 'phylocsf_power_mean', 'Domain', 'is_novel',
            'isoform_switch_q_value', 'gene_q_value', 'IF_overall', 'transcode',
            'leafcutter_sqtl_coloc_match', 'novel_coding_junction_sqtl_coloc_match',
            'riboseq_iso_q', 'riboseq_switch_q',
            'score',
        ])
        for row in reader:
            orf_id = row['ORF_id']
            tx_id  = row['transcript_id']
            TIS_idx = int(row['TIS_pos']) - 1
            TTS_idx = int(row['TTS_pos']) - 1

            loc = orf_locations.get(orf_id)
            location = f"{loc[0]}:{loc[1]}-{loc[2]}" if loc else ''

            gene_id = gene_for_tx.get(tx_id, '')
            gene_name = gene_names.get(gene_id, '')

            phylocsf_bin, power = tx_scores.get(orf_id, ('', ''))
            try:
                phylocsf_bin = 1 if float(phylocsf_bin) > 5 else 0
            except (ValueError, TypeError):
                phylocsf_bin = 0

            sig = collaborator_signatures.get(orf_id)
            # is_novel: the isoform has >=1 splice junction absent from GENCODE.
            # Single-exon isoforms (no junctions) are not novel by this definition.
            orf_junctions = collaborator_junctions.get(orf_id, set())
            is_novel = 1 if any(j not in gencode_junctions for j in orf_junctions) else 0
            transcode = 1 if (sig is not None and sig in study2_signatures) else 0

            domain = 1 if orf_id in pfam_hits else 0
            isoform_switch = switch_q.get(orf_id, 0)
            gene_de = gene_q.get(orf_id, 0)
            leafcutter_match = 1 if gene_name in leafcutter_coloc_genes else 0
            novel_junction_match = 1 if gene_name in novel_junction_coloc_genes else 0
            riboseq_iso_q = classify_riboseq(
                ribo_log2fc.get(orf_id), iso_log2fc.get(orf_id), iso_q_value.get(orf_id)
            )
            riboseq_switch_q = classify_riboseq(
                ribo_log2fc.get(orf_id), dif.get(orf_id), switch_q_value.get(orf_id)
            )
            score = (phylocsf_bin + domain + is_novel + isoform_switch + gene_de
                     + transcode + leafcutter_match + novel_junction_match)

            # ORF_type_ORFanage: canonical CDS from orfanage GTF
            canonical_orfanage = canonical_positions.get(tx_id)
            orf_type_orfanage = classify_orf_type(
                TIS_idx, TTS_idx, canonical_orfanage, tx_id, gene_for_tx, gene_biotypes
            )

            # ORF_type_RiboTIE: GENCODE CDS projected onto the Iso-Seq transcript
            gencode_tx = isoform_to_gencode_tx.get(tx_id, 'novel')
            if (gencode_tx != 'novel'
                    and gencode_tx in gencode_cds_by_tx
                    and tx_id in iso_exons):
                canonical_gencode = project_gencode_cds(
                    iso_exons[tx_id], gencode_cds_by_tx[gencode_tx]
                )
            else:
                canonical_gencode = None
            orf_type_gencode = classify_orf_type(
                TIS_idx, TTS_idx, canonical_gencode, tx_id, gene_for_tx, gene_biotypes
            )

            writer.writerow([
                orf_id,
                gene_name,
                location,
                tx_id,
                row['transcript_len'],
                row['start_codon'],
                row['stop_codon'],
                row['strand'],
                orf_type_orfanage,
                orf_type_gencode,
                phylocsf_bin,
                power,
                domain,
                is_novel,
                isoform_switch,
                gene_de,
                if_overall.get(orf_id, 0),
                transcode,
                leafcutter_match,
                novel_junction_match,
                riboseq_iso_q,
                riboseq_switch_q,
                score,
            ])


if __name__ == '__main__':
    main()
