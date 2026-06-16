#!/usr/bin/env python
"""Trim stop codon (3 nt) from the terminal CDS exon in a RiboTISH GTF.

RiboTISH/ORFannotate includes the stop codon in CDS coordinates, while GENCODE
does not. This script normalizes the ribotish GTF to the GENCODE convention by
trimming 3 nt from the terminal CDS exon of each transcript.

  + strand: subtract 3 from `end` of the CDS exon with the highest end.
  - strand: add 3 to `start` of the CDS exon with the lowest start.

Rows where the trimmed exon would be <= 0 nt are dropped (the terminal CDS
exon was entirely composed of the stop codon).
"""
import argparse
import polars as pl


def read_gtf(path: str, attributes: list[str]) -> pl.DataFrame:
    return (
        pl.read_csv(
            path,
            separator="\t",
            comment_prefix="#",
            schema_overrides={"seqname": pl.String},
            has_header=False,
            new_columns=["seqname", "source", "feature", "start", "end",
                         "score", "strand", "frame", "attributes"],
        )
        .with_columns([
            pl.col("attributes").str.extract(rf'{attr} "([^;]*)";').alias(attr)
            for attr in attributes
        ])
    )


def write_gtf(df: pl.DataFrame, path: str) -> None:
    df.select(["seqname", "source", "feature", "start", "end",
               "score", "strand", "frame", "attributes"]
    ).write_csv(path, separator="\t", include_header=False, quote_style="never")


def trim_stop_codon(gtf: pl.DataFrame) -> pl.DataFrame:
    cds = gtf.filter(pl.col("feature") == "CDS")
    other = gtf.filter(pl.col("feature") != "CDS")

    cds = (
        cds
        .with_columns(
            _max_end=pl.col("end").max().over("transcript_id"),
            _min_start=pl.col("start").min().over("transcript_id"),
        )
        .with_columns(
            _is_terminal=(
                pl.when(pl.col("strand") == "+")
                .then(pl.col("end") == pl.col("_max_end"))
                .otherwise(pl.col("start") == pl.col("_min_start"))
            )
        )
        .with_columns(
            end=pl.when(pl.col("_is_terminal") & (pl.col("strand") == "+"))
                .then(pl.col("end") - 3)
                .otherwise(pl.col("end")),
            start=pl.when(pl.col("_is_terminal") & (pl.col("strand") == "-"))
                .then(pl.col("start") + 3)
                .otherwise(pl.col("start")),
        )
        .drop(["_max_end", "_min_start", "_is_terminal"])
        .filter(pl.col("end") >= pl.col("start"))  # drop exons that were purely stop codon
    )

    return pl.concat([other, cds], how="diagonal").sort(["seqname", "start"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input_gtf", help="RiboTISH GTF (stop codon included in CDS)")
    parser.add_argument("--output_gtf", help="Output GTF (stop codon trimmed from CDS)")
    args = parser.parse_args()

    gtf = read_gtf(args.input_gtf, attributes=["gene_id", "transcript_id"])

    n_before = gtf.filter(pl.col("feature") == "CDS").height
    trimmed = trim_stop_codon(gtf)
    n_after = trimmed.filter(pl.col("feature") == "CDS").height

    dropped = n_before - n_after
    if dropped:
        print(f"Dropped {dropped} CDS rows that were entirely stop codon")
    print(f"CDS rows: {n_before} → {n_after}")

    write_gtf(trimmed, args.output_gtf)
    print(f"Written to {args.output_gtf}")


if __name__ == "__main__":
    main()