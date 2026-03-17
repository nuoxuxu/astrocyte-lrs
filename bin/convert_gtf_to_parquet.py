#!/usr/bin/env python3
import polars as pl
import argparse
import re

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                ).drop("attributes")

def _parse_gtf_attributes(attrs_col: pl.Series) -> pl.DataFrame:
    all_keys = set()
    parsed = []
    for attr_str in attrs_col:
        pairs = re.findall(r'(\w+)\s+"([^"]*)"', attr_str if attr_str else "")
        d = dict(pairs)
        all_keys.update(d.keys())
        parsed.append(d)
    sorted_keys = sorted(all_keys)
    return pl.DataFrame([
        {k: row.get(k) for k in sorted_keys}
        for row in parsed
    ])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_gtf", action = "store", dest="input_gtf")
    parser.add_argument("--output_gtf", action = "store", dest="output_gtf")
    params = parser.parse_args()

    input_gtf = read_gtf(params.input_gtf)
    _attr_df = _parse_gtf_attributes(input_gtf["attributes"])
    _existing_cols = set(input_gtf.columns) - {"attributes"}
    _new_cols = [c for c in _attr_df.columns if c not in _existing_cols]
    input_gtf_parsed = pl.concat(
        [input_gtf.drop("attributes"), _attr_df.select(_new_cols)],
        how="horizontal",
    )

    input_gtf_parsed.write_parquet(params.output_gtf)

if __name__ == "__main__":
    main()
