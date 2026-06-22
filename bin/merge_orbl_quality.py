#!/usr/bin/env python3
"""Merge ORBLv scores, ribotie score, GENCODE label, and classification data into quality metrics TSV.

Usage:
    merge_orbl_quality.py  orf_type_gencode.tsv  orbl_output.tsv
        [--ribotie-csv ribotie_res_merged.csv]
        [--gencode-ribotie gencode/merged/ribotie_res_merged.csv]
        [--classification final_classification.parquet]
        -o quality_metrics.tsv

Joins on ORF_id, appending ORBLv, ribotie_score, in_gencode_ribotie, structural_category,
and gene_name columns. ORFs absent from optional data sources receive 'NA' or False.
"""
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description='Merge ORBLv scores and RiboTIE data into quality TSV.')
    parser.add_argument('orf_type_gencode', help='orf_type_gencode.tsv (ORF_id, transcript_id, ORF_type_GENCODE, ...)')
    parser.add_argument('orbl_output', help='orbl.py output TSV (with --header --extraFields ORF_id,...)')
    parser.add_argument('--ribotie-csv', help='no_minlen RiboTIE merged CSV (for ribotie_score and gene_name)')
    parser.add_argument('--gencode-ribotie', help='GENCODE RiboTIE merged CSV (for in_gencode_ribotie label)')
    parser.add_argument('--classification', help='SQANTI classification parquet file (for structural_category)')
    parser.add_argument('-o', '--output', default='quality_metrics.tsv')
    args = parser.parse_args()

    # Build ORF_id -> ORBLv from orbl output (header present)
    orblv = {}
    with open(args.orbl_output) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            orf_id = row.get('ORF_id', '').strip()
            if orf_id:
                orblv[orf_id] = row.get('ORBLv', 'NA')

    # Build ORF_id -> ribotie_score from no_minlen RiboTIE CSV
    ribotie_score = {}
    if args.ribotie_csv:
        with open(args.ribotie_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                orf_id = row.get('ORF_id', '').strip()
                if orf_id:
                    ribotie_score[orf_id] = row.get('ribotie_score', 'NA')

    # Build protein_seq set from GENCODE RiboTIE CSV for GENCODE matching
    gencode_proteins = set()
    if args.gencode_ribotie:
        with open(args.gencode_ribotie) as f:
            reader = csv.DictReader(f)
            for row in reader:
                protein_seq = row.get('protein_seq', '').strip()
                if protein_seq:
                    gencode_proteins.add(protein_seq)

    # Build ORF_id -> protein_seq from no_minlen RiboTIE CSV
    pacbio_proteins = {}
    # Build ORF_id -> gene_name from no_minlen RiboTIE CSV
    gene_name_map = {}
    if args.ribotie_csv:
        with open(args.ribotie_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                orf_id = row.get('ORF_id', '').strip()
                protein_seq = row.get('protein_seq', '').strip()
                gene_name = row.get('gene_name', '').strip()
                if orf_id and protein_seq:
                    pacbio_proteins[orf_id] = protein_seq
                if orf_id and gene_name:
                    gene_name_map[orf_id] = gene_name

    # Build transcript_id -> structural_category from classification parquet
    structural_category_map = {}
    if args.classification:
        try:
            import polars as pl
            df_class = pl.read_parquet(args.classification)
            for row in df_class.iter_rows(named=True):
                isoform = row.get('isoform', '').strip()
                structural_category = row.get('structural_category', '').strip()
                if isoform and structural_category:
                    structural_category_map[isoform] = structural_category
        except ImportError:
            print("Warning: polars not available, skipping structural_category", file=__import__('sys').stderr)

    with open(args.orf_type_gencode) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')

        # Write header with all columns
        header = ['ORF_id', 'transcript_id', 'ORF_type_GENCODE', 'ORF_type_ORFanage', 'ORBLv']
        if args.ribotie_csv:
            header.append('ribotie_score')
        if args.gencode_ribotie:
            header.append('in_gencode_ribotie')
        if args.classification:
            header.append('structural_category')
        if args.ribotie_csv:
            header.append('gene_name')
        writer.writerow(header)

        for row in reader:
            orf_id = row['ORF_id']
            transcript_id = row['transcript_id']
            out_row = [
                orf_id,
                transcript_id,
                row['ORF_type_GENCODE'],
                row.get('ORF_type_ORFanage', 'NA'),
                orblv.get(orf_id, 'NA'),
            ]

            if args.ribotie_csv:
                out_row.append(ribotie_score.get(orf_id, 'NA'))

            if args.gencode_ribotie:
                protein_seq = pacbio_proteins.get(orf_id, '')
                in_gencode = protein_seq in gencode_proteins if protein_seq else False
                out_row.append(str(in_gencode))

            if args.classification:
                structural_category = structural_category_map.get(transcript_id, 'NA')
                out_row.append(structural_category)

            if args.ribotie_csv:
                out_row.append(gene_name_map.get(orf_id, 'NA'))

            writer.writerow(out_row)


if __name__ == '__main__':
    main()
