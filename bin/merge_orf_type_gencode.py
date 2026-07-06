#!/scratch/nxu/astrocyte-lrs/env/bin/python3
"""Merge ORF_type_ORFanage and ORF_type_RiboTIE labels into the final quality-metrics TSV.

Joins the riboseq quality-metrics file with the orf_type_gencode TSV on ORF_id
and inserts both ORF type columns immediately before ORBLv.
"""
import argparse

import polars as pl

TYPE_COLS = ['ORF_type_ORFanage', 'ORF_type_RiboTIE']


def main():
    parser = argparse.ArgumentParser(description='Merge ORF type labels into quality metrics.')
    parser.add_argument('quality_metrics', help='Riboseq quality metrics TSV (add_riboseq_coverage output)')
    parser.add_argument('orf_type_gencode', help='ORF type TSV (label_orf_type_gencode output)')
    parser.add_argument('-o', '--output', default='final_quality_metrics.tsv')
    args = parser.parse_args()

    metrics = pl.read_csv(args.quality_metrics, separator='\t', infer_schema_length=0)
    labels = pl.read_csv(args.orf_type_gencode, separator='\t').select(['ORF_id'] + TYPE_COLS)

    merged = metrics.join(labels, on='ORF_id', how='left')

    cols = [c for c in merged.columns if c not in TYPE_COLS]
    insert_pos = cols.index('ORBLv')
    ordered = cols[:insert_pos] + TYPE_COLS + cols[insert_pos:]
    merged.select(ordered).write_csv(args.output, separator='\t')


if __name__ == '__main__':
    main()
