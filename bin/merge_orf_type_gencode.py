#!/scratch/nxu/astrocyte-lrs/env/bin/python3
"""Merge ORF_type_GENCODE labels into the final quality-metrics TSV.

Joins the riboseq quality-metrics file with the orf_type_gencode TSV on ORF_id
and inserts the ORF_type_GENCODE column immediately after ORF_type_ORFanage.
"""
import argparse

import polars as pl


def main():
    parser = argparse.ArgumentParser(description='Merge ORF_type_GENCODE into quality metrics.')
    parser.add_argument('quality_metrics', help='Riboseq quality metrics TSV (add_riboseq_coverage output)')
    parser.add_argument('orf_type_gencode', help='ORF type GENCODE TSV (label_orf_type_gencode output)')
    parser.add_argument('-o', '--output', default='final_quality_metrics.tsv')
    args = parser.parse_args()

    metrics = pl.read_csv(args.quality_metrics, separator='\t', infer_schema_length=0)
    labels = pl.read_csv(args.orf_type_gencode, separator='\t').select(['ORF_id', 'ORF_type_GENCODE'])

    merged = metrics.join(labels, on='ORF_id', how='left')

    cols = merged.columns
    insert_pos = cols.index('ORF_type_ORFanage') + 1
    ordered = cols[:insert_pos] + ['ORF_type_GENCODE'] + [c for c in cols[insert_pos:] if c != 'ORF_type_GENCODE']
    merged.select(ordered).write_csv(args.output, separator='\t')


if __name__ == '__main__':
    main()
