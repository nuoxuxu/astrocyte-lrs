#!/bin/env python
import polars as pl
import argparse
from src.utils import read_gtf


def main():
    parser = argparse.ArgumentParser(description='Add lncRNA ORF labels to RiboTIE CSV.')
    parser.add_argument('ribotie_csv', help='RiboTIE predictions CSV')
    parser.add_argument('final_classification', help='SQANTI3 final_classification.parquet')
    parser.add_argument('annotation_gtf', help='GENCODE annotation GTF')
    parser.add_argument('-o', '--output', default='ribotie_cpm1_3sample_with_lncRNA.csv',
                        help='Output CSV file (default: ribotie_cpm1_3sample_with_lncRNA.csv)')
    args = parser.parse_args()

    ribotie_cpm1_3sample = pl.read_csv(args.ribotie_csv)
    classification = pl.read_parquet(args.final_classification)
    gencode_gtf = read_gtf(args.annotation_gtf, attributes=["gene_id", "gene_name", "transcript_id", "transcript_name", "transcript_type", "gene_type"])

    ribotie_cpm1_3sample\
        .join(classification["isoform", "associated_gene"], left_on="transcript_id", right_on="isoform")\
        .join(
            gencode_gtf.filter(pl.col("feature")=="gene").select("gene_id", "gene_type"), left_on="associated_gene", right_on="gene_id", how="left"
        )\
        .with_columns(
            pl.when(pl.col("gene_type")=="lncRNA").then(pl.lit("lncRNA ORF")).otherwise(pl.col("ORF_type")).alias("ORF_type")
        )\
        .drop("associated_gene", "gene_type")\
        .write_csv(args.output)


if __name__ == '__main__':
    main()
