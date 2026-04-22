#!/bin/env python
import argparse
import polars as pl
from src.utils import read_gtf

parser = argparse.ArgumentParser(description="Filter RiboTIE GTF results by CPM threshold CSV")
parser.add_argument("--input_gtf", required=True, help="Input RiboTIE GTF file")
parser.add_argument("--ribotie_cpm1_3sample", required=True, help="CSV with ORF_id column listing ORFs passing CPM threshold")
parser.add_argument("--output_gtf", required=True, help="Output filtered GTF file")
args = parser.parse_args()

ribotie_gtf = read_gtf(args.input_gtf, attributes=["transcript_id", "ORF_id"])
ribotie_cpm1_3sample = pl.read_csv(args.ribotie_cpm1_3sample)

non_coding = ribotie_gtf.filter(pl.col("ORF_id").is_null())
coding = ribotie_gtf.filter(pl.col("transcript_id").is_in(ribotie_cpm1_3sample.get_column("ORF_id").unique().to_list()))
filtered = pl.concat([non_coding, coding])

filtered\
    .drop(["transcript_id", "ORF_id"])\
    .write_csv(args.output_gtf, separator="\t", include_header=False, quote_style="never")
