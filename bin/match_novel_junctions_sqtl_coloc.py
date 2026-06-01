#!/usr/bin/env python3
"""
Cross-reference novel coding junction sQTLs with BigBrain pre-computed COLOC results.

Joins novel_coding_junction_sqtl_matches.tsv with BigBrain COLOC results (PP.H4 >= 0.5)
on the shared feature column (chrN_start_end format). The IF column (max IF_overall
across transcripts carrying each junction) is carried through from the input.
Output has one row per junction x disease x GWAS locus combination.
"""

import argparse

import polars as pl


def main():
    parser = argparse.ArgumentParser(
        description="Cross-reference novel coding junction sQTLs with BigBrain COLOC results."
    )
    parser.add_argument("sqtl_matches", help="Novel coding junction sQTL matches TSV")
    parser.add_argument("coloc_file", help="BigBrain COLOC results TSV (PP.H4 >= 0.5)")
    parser.add_argument("-o", "--output", default="novel_coding_junction_sqtl_coloc_matches.tsv",
                        help="Output TSV file (default: novel_coding_junction_sqtl_coloc_matches.tsv)")
    args = parser.parse_args()

    sqtl = pl.read_csv(args.sqtl_matches, separator="\t", null_values="NA")
    coloc = pl.read_csv(args.coloc_file, separator="\t", null_values="NA")
    result = sqtl.join(coloc, on="feature", how="inner")

    result.write_csv(args.output, separator="\t")
    print(f"Matched rows:      {len(result)}")
    print(f"Unique features:   {result['feature'].n_unique()}")
    print(f"Diseases present:  {sorted(result['disease'].unique().to_list())}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
