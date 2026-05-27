#!/usr/bin/env python3

import argparse
import polars as pl
from src.utils import read_fasta, write_fasta


def read_gtf_transcript_ids(gtf_file):
    return (
        pl.read_csv(
            gtf_file,
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            has_header=False,
            new_columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"],
            schema_overrides={"seqname": pl.String, "column_8": pl.String},
        )
        .filter(pl.col("feature") == "transcript")
        .with_columns(
            pl.col("attributes").str.extract(r'transcript_id "([^;]*)";').alias("transcript_id")
        )
        .with_columns(
            pl.col("transcript_id").map_elements(
                lambda x: x.rsplit("_", 1)[0] if "_" in x else x,
                return_dtype=pl.String,
            ).alias("base_id")
        )
        .select(["transcript_id", "base_id"])
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("supplemented_gtf")
    parser.add_argument("input_fasta")
    parser.add_argument("input_expression")
    parser.add_argument("--output_fasta", default="supplemented_collaborator.fasta")
    parser.add_argument("--output_expression", default="supplemented_collaborator_expression.parquet")
    args = parser.parse_args()

    tx_df = read_gtf_transcript_ids(args.supplemented_gtf)

    expression = pl.read_parquet(args.input_expression)
    (
        tx_df
        .join(expression, left_on="base_id", right_on="tname", how="left")
        .drop("base_id")
        .rename({"transcript_id": "tname"})
        .write_parquet(args.output_expression)
    )

    fasta_df = read_fasta(args.input_fasta)
    result_fasta = (
        tx_df
        .join(fasta_df, left_on="base_id", right_on="transcript_id", how="left")
        .drop("base_id")
    )
    write_fasta(result_fasta, args.output_fasta)


if __name__ == "__main__":
    main()
