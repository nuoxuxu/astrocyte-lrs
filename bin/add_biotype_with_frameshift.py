#!/usr/bin/env python3
"""Add biotype_with_frameshift column to quality_metrics.tsv for ORBLq calculation.

For each ORF, reprojects the GENCODE canonical CDS onto the Iso-Seq transcript
coordinate system (same logic as label_orf_type_gencode.py) and computes the
frame offset between the ORF TIS and the GENCODE canonical TIS:

    frame = (TIS_idx - canonical_TIS_idx) % 3

Overlapping types (uoORF, intORF, doORF) are annotated as e.g. uoORF+1 or uoORF+2;
frame 0 (in-frame with canonical CDS) yields the bare type name without a suffix.
Types with no frameshift variant (uORF, dORF, lncRNA-ORF) are passed through directly.
All other types (annotated CDS, N/C-terminal extension/truncation, varRNA-ORF) get NA.

BiotypesWithFS = ['uORF', 'uoORF+1', 'uoORF+2', 'intORF+1', 'intORF+2',
                  'doORF+1', 'doORF+2', 'dORF', 'lncRNA-ORF']

Usage:
    add_biotype_with_frameshift.py quality_metrics.tsv ribotie_csv annotation_gtf \
        orfanage_gtf final_classification -o output.tsv
"""
import argparse
import csv
import re

import polars as pl

DIRECT = {'uORF': 'uORF', 'dORF': 'dORF', 'lncRNA-ORF': 'lncRNA-ORF'}
NEEDS_FRAME = {'uoORF', 'intORF', 'doORF'}


def parse_gtf_attributes(attr_str):
    attrs = {}
    for match in re.finditer(r'(\w+) "([^"]*)"', attr_str):
        attrs[match.group(1)] = match.group(2)
    return attrs


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


def build_frame_map(ribotie_csv, iso_exons, gencode_cds_by_tx,
                    isoform_to_gencode_tx, gene_for_tx, gene_to_coding_txs):
    """Return ORF_id -> frame (int 0/1/2 or None) for all overlapping ORF types."""
    frame_map = {}
    with open(ribotie_csv) as f:
        for row in csv.DictReader(f):
            orf_id = row['ORF_id']
            tx_id = row['transcript_id']
            tis_idx = int(row['TIS_pos']) - 1

            # Collect all GENCODE coding transcripts to consider (same logic as
            # label_orf_type_gencode.py): SQANTI-assigned transcript + all transcripts
            # of each gene in associated_gene, splitting compound fusion gene IDs.
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
                        dist = abs(c[0] - tis_idx)
                        if dist < best_dist:
                            best, best_dist = c, dist
                canonical = best

            frame_map[orf_id] = (tis_idx - canonical[0]) % 3 if canonical is not None else None
    return frame_map


def main():
    parser = argparse.ArgumentParser(
        description='Add biotype_with_frameshift to quality_metrics.tsv.'
    )
    parser.add_argument('quality_metrics')
    parser.add_argument('ribotie_csv')
    parser.add_argument('annotation_gtf')
    parser.add_argument('orfanage_gtf')
    parser.add_argument('final_classification')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    _gene_biotypes, gencode_cds_by_tx, _gene_names, gene_to_coding_txs = load_gencode_data(args.annotation_gtf)
    iso_exons = load_orfanage_iso_exons(args.orfanage_gtf)

    clf = pl.read_parquet(args.final_classification).select(
        ['isoform', 'associated_gene', 'associated_transcript']
    )
    isoform_to_gencode_tx = dict(zip(clf['isoform'], clf['associated_transcript']))
    gene_for_tx = dict(zip(clf['isoform'], clf['associated_gene']))

    frame_map = build_frame_map(
        args.ribotie_csv, iso_exons, gencode_cds_by_tx,
        isoform_to_gencode_tx, gene_for_tx, gene_to_coding_txs
    )

    with open(args.quality_metrics) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow([
            'ORF_id', 'transcript_id', 'ORF_type_GENCODE', 'ORBLv',
            'biotype_with_frameshift',
        ])
        for row in reader:
            orf_id = row['ORF_id']
            orf_type = row['ORF_type_GENCODE']

            if orf_type in DIRECT:
                biotype_fs = DIRECT[orf_type]
            elif orf_type in NEEDS_FRAME:
                frame = frame_map.get(orf_id)
                biotype_fs = f'{orf_type}+{frame}' if frame in (1, 2) else orf_type
            else:
                biotype_fs = 'NA'

            writer.writerow([
                orf_id,
                row['transcript_id'],
                orf_type,
                row['ORBLv'],
                biotype_fs,
            ])


if __name__ == '__main__':
    main()
