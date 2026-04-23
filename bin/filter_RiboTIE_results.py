#!/bin/env python
import argparse
import polars as pl
from src.utils import read_gtf, read_fasta, write_fasta

def main():
    parser = argparse.ArgumentParser(description="Filter RiboTIE GTF results by CPM threshold CSV")
    parser.add_argument("--input_gtf", required=True, help="Input RiboTIE GTF file")
    parser.add_argument("--ribotie_cpm1_3sample", required=True, help="CSV with ORF_id column listing ORFs passing CPM threshold")
    parser.add_argument("--input_fasta", required=True, help="Input FASTA file of final transcripts")
    parser.add_argument("--input_expression", required=True, help="Input Parquet file of final expression")
    parser.add_argument("--output_gtf", required=True, help="Output filtered GTF file")
    parser.add_argument("--output_fasta", required=True, help="Output filtered FASTA file")
    parser.add_argument("--output_expression", required=True, help="Output filtered expression Parquet file")
    args = parser.parse_args()

    ribotie_gtf = read_gtf(args.input_gtf, attributes=["transcript_id", "ORF_id"])
    ribotie_cpm1_3sample = pl.read_csv(args.ribotie_cpm1_3sample)

    non_coding = ribotie_gtf.filter(pl.col("ORF_id").is_null())
    coding = ribotie_gtf.filter(pl.col("transcript_id").is_in(ribotie_cpm1_3sample.get_column("ORF_id").unique().to_list()))
    filtered = pl.concat([non_coding, coding])

    filtered\
        .drop(["transcript_id", "ORF_id"])\
        .write_csv(args.output_gtf, separator="\t", include_header=False, quote_style="never")

    tx_list = filtered\
        .filter(pl.col("feature") == "transcript")\
        .with_columns(
            pl.col("transcript_id").map_elements(lambda x: x.split("_")[0], return_dtype=pl.String).alias("transcript_id_base")
        )\
        .select(["transcript_id", "transcript_id_base"])

    final_fasta = read_fasta(args.input_fasta)
    final_expression = pl.read_parquet(args.input_expression)

    tx_list\
        .join(final_expression, left_on="transcript_id_base", right_on="tname", how="left")\
        .drop("transcript_id_base")\
        .rename({"transcript_id": "tname"})\
        .write_parquet(args.output_expression)

    ribotie_filtered_final_fasta = tx_list\
        .join(final_fasta, left_on="transcript_id_base", right_on="transcript_id", how="left")\
        .drop("transcript_id_base")

    write_fasta(ribotie_filtered_final_fasta, args.output_fasta)

if __name__ == "__main__":
    main()
