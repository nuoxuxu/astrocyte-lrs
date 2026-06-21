#!/usr/bin/env python3
"""Merge ORBLv scores from orbl.py output into the orf_type_gencode TSV.

Usage:
    merge_orbl_quality.py  orf_type_gencode.tsv  orbl_output.tsv  -o quality_metrics.tsv

Joins on ORF_id, appending the ORBLv column. ORFs absent from the ORBL output
receive 'NA'.
"""
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description='Merge ORBLv scores into quality TSV.')
    parser.add_argument('orf_type_gencode', help='orf_type_gencode.tsv (ORF_id, transcript_id, ORF_type_GENCODE)')
    parser.add_argument('orbl_output', help='orbl.py output TSV (with --header --extraFields ORF_id,...)')
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

    with open(args.orf_type_gencode) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['ORF_id', 'transcript_id', 'ORF_type_GENCODE', 'ORBLv'])
        for row in reader:
            orf_id = row['ORF_id']
            writer.writerow([
                orf_id,
                row['transcript_id'],
                row['ORF_type_GENCODE'],
                orblv.get(orf_id, 'NA'),
            ])


if __name__ == '__main__':
    main()
