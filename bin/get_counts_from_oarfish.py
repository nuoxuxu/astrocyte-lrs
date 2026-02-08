#!/bin/env python
import polars as pl
from pathlib import Path
import polars.selectors as cs
from src.utils import read_gtf
import argparse

parser = argparse.ArgumentParser(description="Get counts from oarfish output")
parser.add_argument("--min_reads", type=int, default=5, help="Minimum number of reads to consider a transcript as expressed")
parser.add_argument("--min_n_sample", type=int, default=2, help="Minimum number of samples to consider a transcript as expressed")
params = parser.parse_args()

def get_exp_tx(expression, min_reads=5, min_n_sample=2):
    return expression\
        .with_columns(
            cs.numeric() > min_reads
        )\
        .filter(
            pl.sum_horizontal(cs.boolean()) > min_n_sample
        )\
        ["tname"].to_list()

def main():
    tx_count = pl.scan_csv(list(Path(".").glob("*.quant")), include_file_paths="path", separator="\t")

    tx_count = tx_count\
        .collect()\
        .with_columns(
            pl.col("path").str.replace("_5p--IsoSeqX_3p.quant", "").map_elements(lambda x: x.split("_")[1], return_dtype=pl.String)
        )\
        .with_columns(pl.col("path").map_elements(lambda x: x.upper(), return_dtype=pl.String))\
        .rename({"path": "sample"})

    # primer_to_sample = pl.read_csv("data/primer_to_sample.csv")

    # tx_count\
    #     .join(
    #         primer_to_sample, 
    #         left_on="sample", right_on="Index primer", 
    #         how="left"
    #     )

    tx_count = tx_count.pivot("sample", index="tname", values="num_reads")

    isoforms_to_keep = get_exp_tx(tx_count, min_reads=params.min_reads, min_n_sample=params.min_n_sample)

    outputDir = f"final_expression_{params.min_reads}_{params.min_n_sample}"

    tx_count.filter(pl.col("tname").is_in(isoforms_to_keep)).write_parquet(f"{outputDir}/final_expression.parquet")

    classification = pl.read_csv("default_RulesFilter_result_classification.txt", separator="\t")\
        .filter(
            pl.col("filter_result")=="Isoform",
            pl.col("isoform").is_in(isoforms_to_keep)
        )

    classification.drop("filter_result").write_parquet(f"{outputDir}/final_classification.parquet")

    read_gtf("default.filtered.gtf")\
        .filter(pl.col("transcript_id").is_in(isoforms_to_keep))\
        .drop("transcript_id")\
        .write_csv(f"{outputDir}/final_transcripts.gtf", separator="\t", include_header=False, quote_style="never")
    
if __name__ == "__main__":
    main()