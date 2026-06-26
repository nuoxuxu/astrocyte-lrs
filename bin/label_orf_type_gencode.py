#!/scratch/nxu/astrocyte-lrs/env/bin/python3
"""Assign ORF_type_GENCODE labels to RiboTIE ORF predictions.

For each ORF in the RiboTIE CSV, finds the canonical GENCODE CDS via a
three-tier hierarchy:
  1. orfanage_template attribute on the matching transcript in orfanage.gtf
  2. associated_transcript from SQANTI3 final_classification
  3. MANE_Select transcript for the associated_gene

Outputs a TSV with columns: ORF_id, transcript_id, ORF_type_ORFanage, ORF_type_GENCODE.
"""
import argparse
import csv
import re

import polars as pl


def parse_gtf_attributes(attr_str):
    attrs = {}
    for match in re.finditer(r'(\w+) "([^"]*)"', attr_str):
        attrs[match.group(1)] = match.group(2)
    return attrs


def get_attr(attr_str, key):
    match = re.search(re.escape(key) + r' "([^"]*)"', attr_str)
    return match.group(1) if match else None


def strip_version(enst_id):
    """Strip ENSEMBL version suffix: ENST00000361390.2 -> ENST00000361390."""
    return enst_id.split('.')[0] if enst_id else enst_id


_FUSION_RE = re.compile(r'^ENSG\d+\.\d+_ENSG\d+\.\d+$')


def load_gencode_data(gtf_path):
    """Parse GENCODE GTF; return (gene_biotypes, cds_by_tx, gene_names, gene_to_mane_select).

    cds_by_tx is keyed by version-stripped transcript IDs.
    gene_to_mane_select maps stripped gene_id -> stripped transcript_id for MANE_Select transcripts.
    """
    tx_cds = {}
    gene_biotypes = {}
    gene_names = {}
    gene_to_mane_select = {}
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
                if g and attrs.get('gene_type'):
                    gene_biotypes[g] = attrs['gene_type']
                if g and attrs.get('gene_name'):
                    gene_names[g] = attrs['gene_name']
            elif feat == 'transcript':
                attrs = parse_gtf_attributes(fields[8])
                tx_id = attrs.get('transcript_id')
                gene_id = attrs.get('gene_id')
                if tx_id and gene_id and 'MANE_Select' in fields[8]:
                    gene_to_mane_select[strip_version(gene_id)] = strip_version(tx_id)
            elif feat == 'CDS':
                chrom, strand = fields[0], fields[6]
                start, end = int(fields[3]), int(fields[4])
                attrs = parse_gtf_attributes(fields[8])
                tx_id = attrs.get('transcript_id')
                if not tx_id:
                    continue
                key = strip_version(tx_id)
                if key not in tx_cds:
                    tx_cds[key] = {'chrom': chrom, 'strand': strand, 'segs': []}
                tx_cds[key]['segs'].append((start, end))
    cds_by_tx = {tx: {'strand': v['strand'], 'segs': v['segs']} for tx, v in tx_cds.items()}
    return gene_biotypes, cds_by_tx, gene_names, gene_to_mane_select


def load_orfanage_iso_exons(gtf_path):
    """Parse orfanage numbered-exons GTF; return iso_exons dict: transcript_id -> {strand, exons}."""
    exons = {}
    strands = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            strand = fields[6]
            start, end = int(fields[3]), int(fields[4])
            attrs = parse_gtf_attributes(fields[8])
            tx_id = attrs.get('transcript_id')
            if not tx_id:
                continue
            strands[tx_id] = strand
            exons.setdefault(tx_id, []).append((start, end))
    return {tx: {'strand': strands[tx], 'exons': exons[tx]} for tx in exons}


def load_orfanage_templates(gtf_path):
    """Parse plain orfanage.gtf for transcript-level orfanage_template attributes.

    Returns {transcript_id: stripped_enst_id}.
    """
    templates = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            attrs = parse_gtf_attributes(fields[8])
            tx_id = attrs.get('transcript_id')
            template = attrs.get('orfanage_template')
            if tx_id and template:
                templates[tx_id] = strip_version(template)
    return templates


def project_gencode_cds(iso_exon_data, cds_data):
    """Project a GENCODE CDS onto an Iso-Seq transcript's coordinate system.

    Returns (TIS_idx, LTS_idx, TTS_idx) or None if the GENCODE TIS is not
    within any Iso-Seq exon.
    """
    strand = iso_exon_data['strand']
    exons = iso_exon_data['exons']
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
            tis_idx = cum + (tis - es)
            break
        elif strand == '-' and es <= tis <= ee:
            tis_idx = cum + (ee - tis)
            break
        cum += ee - es + 1

    if tis_idx == -1:
        return None

    cds_len = 0
    for cs, ce in cds_segs:
        for es, ee in exons:
            ovl_s, ovl_e = max(cs, es), min(ce, ee)
            if ovl_s <= ovl_e:
                cds_len += ovl_e - ovl_s + 1

    tts_idx = tis_idx + cds_len
    return (tis_idx, tts_idx - 1, tts_idx)


def classify_orf_type(TIS_idx, TTS_idx, canonical, tx_id, gene_for_tx, gene_biotypes):
    """Classify ORF type using the transcript_transformer decision tree."""
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


def resolve_canonical_enst(tx_id, orfanage_templates, assoc_tx_raw, gene_id_raw,
                            gencode_cds_by_tx, gene_to_mane_select):
    """Return (canonical_enst_stripped, source) or (no_template_label, 'special').

    source is one of: 'orfanage_template', 'assoc_tx', 'mane_select', 'special'.
    When source == 'special', the returned string is the final ORF_type_GENCODE label.
    """
    # Tier 1: orfanage_template
    template = orfanage_templates.get(tx_id)
    if template:
        return template, 'orfanage_template'

    # Tier 2: SQANTI3 associated_transcript
    if assoc_tx_raw and assoc_tx_raw != 'novel':
        return strip_version(assoc_tx_raw), 'assoc_tx'

    # Tier 3: MANE_Select for the associated gene
    if not gene_id_raw:
        return 'no_template_unknown_gene', 'special'

    # Fusion: two ENSG IDs joined by underscore
    if _FUSION_RE.match(gene_id_raw):
        return 'no_template_fusion', 'special'

    # Novel gene: not in GENCODE at all
    if gene_id_raw.startswith('novelGene_'):
        return 'no_template_novel_gene', 'special'

    # Generic underscore pattern not matching the above → treat as fusion per design
    if '_' in gene_id_raw:
        return 'no_template_fusion', 'special'

    gene_stripped = strip_version(gene_id_raw)
    mane_tx = gene_to_mane_select.get(gene_stripped)
    if mane_tx is None:
        return 'no_template_no_MANE_Select', 'special'

    return mane_tx, 'mane_select'


def main():
    parser = argparse.ArgumentParser(description='Add ORF_type_GENCODE labels to RiboTIE CSV.')
    parser.add_argument('ribotie_csv', help='RiboTIE predictions CSV')
    parser.add_argument('annotation_gtf', help='GENCODE annotation GTF')
    parser.add_argument('orfanage_gtf', help='Orfanage numbered-exons GTF (orfanage_numbered_exons.gtf)')
    parser.add_argument('final_classification', help='SQANTI3 final_classification.parquet')
    parser.add_argument('orfanage_plain_gtf', help='Plain orfanage.gtf with orfanage_template attributes')
    parser.add_argument('-o', '--output', default='orf_type_gencode.tsv')
    args = parser.parse_args()

    gene_biotypes, gencode_cds_by_tx, _gene_names, gene_to_mane_select = load_gencode_data(args.annotation_gtf)
    iso_exons = load_orfanage_iso_exons(args.orfanage_gtf)
    orfanage_templates = load_orfanage_templates(args.orfanage_plain_gtf)

    clf = pl.read_parquet(args.final_classification).select(
        ['isoform', 'associated_gene', 'associated_transcript']
    )
    gene_for_tx = dict(zip(clf['isoform'], clf['associated_gene']))
    isoform_to_gencode_tx = dict(zip(clf['isoform'], clf['associated_transcript']))

    with open(args.ribotie_csv) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in)
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['ORF_id', 'transcript_id', 'ORF_type_ORFanage', 'ORF_type_GENCODE'])
        for row in reader:
            orf_id = row['ORF_id']
            tx_id = row['transcript_id']
            orf_type_orfanage = row.get('ORF_type', 'NA')
            TIS_idx = int(row['TIS_pos']) - 1
            TTS_idx = int(row['TTS_pos']) - 1

            assoc_tx_raw = isoform_to_gencode_tx.get(tx_id)
            gene_id_raw = gene_for_tx.get(tx_id)

            # tx_id not in classification file at all
            if assoc_tx_raw is None and gene_id_raw is None:
                writer.writerow([orf_id, tx_id, orf_type_orfanage, 'no_template_unknown_gene'])
                continue

            canonical_enst, source = resolve_canonical_enst(
                tx_id, orfanage_templates, assoc_tx_raw, gene_id_raw,
                gencode_cds_by_tx, gene_to_mane_select
            )

            if source == 'special':
                writer.writerow([orf_id, tx_id, orf_type_orfanage, canonical_enst])
                continue

            # Resolve canonical CDS; map non-coding ENST to appropriate label
            if canonical_enst not in gencode_cds_by_tx:
                if source == 'orfanage_template':
                    label = 'no_template_orfanage_template_non_coding'
                elif source == 'assoc_tx':
                    label = 'no_template_associated_transcript_non_coding'
                else:  # mane_select
                    label = 'no_template_associated_MANE_Select_non_coding'
                writer.writerow([orf_id, tx_id, orf_type_orfanage, label])
                continue

            canonical = None
            if tx_id in iso_exons:
                canonical = project_gencode_cds(iso_exons[tx_id], gencode_cds_by_tx[canonical_enst])

            orf_type_gencode = classify_orf_type(
                TIS_idx, TTS_idx, canonical, tx_id, gene_for_tx, gene_biotypes
            )
            writer.writerow([orf_id, tx_id, orf_type_orfanage, orf_type_gencode])


if __name__ == '__main__':
    main()
