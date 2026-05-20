#!/bin/env python
import argparse
import types
import polars as pl
from src.utils import read_fasta, write_fasta

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", quote_char=None, schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", quote_char=None, schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                ).drop("attributes")
    
def main():
    parser = argparse.ArgumentParser(description="Filter RiboTIE GTF results by CPM threshold CSV")
    parser.add_argument("--input_gtf", required=True, help="Input RiboTIE GTF file")
    parser.add_argument("--input_fasta", required=True, help="Input FASTA file of final transcripts")
    parser.add_argument("--input_expression", required=True, help="Input Parquet file of final expression")
    parser.add_argument("--input_classification", required=True, help="Input Parquet file of final classification")
    parser.add_argument("--output_fasta", required=True, help="Output filtered FASTA file")
    parser.add_argument("--output_expression", required=True, help="Output filtered expression Parquet file")
    parser.add_argument("--output_classification", required=True, help="Output filtered classification Parquet file")
    args = parser.parse_args()

    ribotie_gtf = read_gtf(args.input_gtf, attributes=["transcript_id"])

    tx_list = ribotie_gtf\
        .filter(pl.col("feature") == "transcript")\
        .with_columns(
            pl.col("transcript_id").map_elements(lambda x: x.split("_")[0], return_dtype=pl.String).alias("transcript_id_base")
        )\
        .select(["transcript_id", "transcript_id_base"])

    final_fasta = read_fasta(args.input_fasta)
    final_expression = pl.read_parquet(args.input_expression)
    final_classification = pl.read_parquet(args.input_classification)

    tx_list\
        .join(final_expression, left_on="transcript_id_base", right_on="tname", how="left")\
        .drop("transcript_id_base")\
        .rename({"transcript_id": "tname"})\
        .write_parquet(args.output_expression)

    ribotie_filtered_final_fasta = tx_list\
        .join(final_fasta, left_on="transcript_id_base", right_on="transcript_id", how="left")\
        .drop("transcript_id_base")

    write_fasta(ribotie_filtered_final_fasta, args.output_fasta)

    final_classification\
        .filter(pl.col("isoform").is_in(tx_list.get_column("transcript_id_base").to_list()))\
        .write_parquet(args.output_classification)

if __name__ == "__main__":
    main()