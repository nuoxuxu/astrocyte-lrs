#!/bin/env python
import json
from pathlib import Path
import numpy as np
import polars as pl
import argparse

def get_num_reads_flnc_poly(path_to_json):
    with open(path_to_json, 'r') as file:
        data = json.load(file)
    return [d for d in data["attributes"] if d.get('id') == "num_reads_flnc_polya"][0]['value']

def get_num_reads_flnc_poly_from_fofn(fofn):
    with open(fofn, 'r') as file:
        flist = [line.strip() for line in file]
    return np.array([get_num_reads_flnc_poly(flist[i]) for i in range(len(flist))]).sum()

def main():
    parser = argparse.ArgumentParser(
        description='Calculate the total number of Full-Length Non-Chimeric Reads with Poly-A Tail reads for each sample based on the FOFN files.'
    )
    parser.add_argument('--path_primer_to_sample', help='Path to the CSV file that maps primers to sample names.')
    parser.add_argument('--output_path', help='Path to the output CSV file that will contain the total number of reads for each sample.')
    args = parser.parse_args()

    arr = np.array([get_num_reads_flnc_poly_from_fofn(fofn) for fofn in list(Path("./").glob("*.fofn"))])

    primer_to_sample = pl.read_csv(args.path_primer_to_sample)

    df = pl.DataFrame(
        {
            "sample": [str(file_name) for file_name in list(Path("./").glob("*.fofn"))],
            "num_reads_flnc_poly": arr
        }
    )

    df\
        .with_columns(pl.col("sample").str.replace("_5p--IsoSeqX_3p.fofn", ""))\
        .with_columns(pl.col("sample").str.replace("IsoSeqX_", ""))\
        .with_columns(pl.col("sample").str.to_uppercase())\
        .join(primer_to_sample, left_on="sample", right_on="Index primer")\
        .write_csv(args.output_path, include_header=True, quote_style="never")
    
if __name__ == "__main__":
    main()