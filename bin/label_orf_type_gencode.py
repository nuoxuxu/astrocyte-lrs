#!/usr/bin/env python3
"""Assign ORF_type_GENCODE labels to RiboTIE ORF predictions.

For each ORF in the RiboTIE CSV, projects the GENCODE canonical CDS onto the
Iso-Seq transcript coordinate system and applies the same classification tree
used in make_summary_table.py.

Outputs a TSV with columns: ORF_id, transcript_id, ORF_type_GENCODE.
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


def load_gencode_data(gtf_path):
    """Parse GENCODE GTF; return (gene_biotypes, cds_by_tx, gene_names, gene_to_coding_txs)."""
    tx_cds = {}
    tx_to_gene = {}
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
                if g and attrs.get('gene_type'):
                    gene_biotypes[g] = attrs['gene_type']
                if g and attrs.get('gene_name'):
                    gene_names[g] = attrs['gene_name']
            elif feat == 'CDS':
                chrom, strand = fields[0], fields[6]
                start, end = int(fields[3]), int(fields[4])
                attrs = parse_gtf_attributes(fields[8])
                tx_id = attrs.get('transcript_id')
                gene_id = attrs.get('gene_id')
                if not tx_id:
                    continue
                if tx_id not in tx_cds:
                    tx_cds[tx_id] = {'chrom': chrom, 'strand': strand, 'segs': []}
                tx_cds[tx_id]['segs'].append((start, end))
                if gene_id and tx_id not in tx_to_gene:
                    tx_to_gene[tx_id] = gene_id
    cds_by_tx = {tx: {'strand': v['strand'], 'segs': v['segs']} for tx, v in tx_cds.items()}
    gene_to_coding_txs = {}
    for tx_id, gene_id in tx_to_gene.items():
        gene_to_coding_txs.setdefault(gene_id, []).append(tx_id)
    return gene_biotypes, cds_by_tx, gene_names, gene_to_coding_txs


def load_orfanage_iso_exons(gtf_path):
    """Parse orfanage GTF; return iso_exons dict: transcript_id -> {strand, exons}."""
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


def main():
    parser = argparse.ArgumentParser(description='Add ORF_type_GENCODE labels to RiboTIE CSV.')
    parser.add_argument('ribotie_csv', help='RiboTIE predictions CSV')
    parser.add_argument('annotation_gtf', help='GENCODE annotation GTF')
    parser.add_argument('orfanage_gtf', help='Orfanage GTF (orfanage_numbered_exons.gtf)')
    parser.add_argument('final_classification', help='SQANTI3 final_classification.parquet')
    parser.add_argument('-o', '--output', default='orf_type_gencode.tsv')
    args = parser.parse_args()

    gene_biotypes, gencode_cds_by_tx, _gene_names, gene_to_coding_txs = load_gencode_data(args.annotation_gtf)
    iso_exons = load_orfanage_iso_exons(args.orfanage_gtf)

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

            # Collect all GENCODE coding transcripts to consider:
            #   1. The SQANTI-assigned transcript (if it has a CDS) — handles fusions
            #      where the compound associated_gene won't match any gene_to_coding_txs key
            #   2. All coding transcripts of each gene in associated_gene (handles
            #      compound/fusion gene IDs by splitting on '_ENSG')
            # Pick the one whose projected TIS is closest to the ORF TIS.
            canonical = None
            if tx_id in iso_exons:
                candidates = set()
                assoc_tx = isoform_to_gencode_tx.get(tx_id, 'novel')
                if assoc_tx != 'novel' and assoc_tx in gencode_cds_by_tx:
                    candidates.add(assoc_tx)
                gene_id = gene_for_tx.get(tx_id, '')
                for part in (gene_id or '').split('_ENSG'):
                    gid = part if part.startswith('ENSG') else ('ENSG' + part if part else '')
                    candidates.update(gene_to_coding_txs.get(gid, []))
                best, best_dist = None, float('inf')
                for alt_tx in candidates:
                    c = project_gencode_cds(iso_exons[tx_id], gencode_cds_by_tx[alt_tx])
                    if c is not None:
                        dist = abs(c[0] - TIS_idx)
                        if dist < best_dist:
                            best, best_dist = c, dist
                canonical = best

            orf_type_gencode = classify_orf_type(
                TIS_idx, TTS_idx, canonical, tx_id, gene_for_tx, gene_biotypes
            )
            writer.writerow([orf_id, tx_id, orf_type_orfanage, orf_type_gencode])


if __name__ == '__main__':
    main()
