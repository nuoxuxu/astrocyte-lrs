#!/usr/bin/env python3
"""
Cross-reference novel coding junction sQTLs with BigBrain pre-computed COLOC results.

Joins novel_coding_junction_sqtl_matches.tsv with BigBrain COLOC results (PP.H4 >= 0.5)
on the shared feature column (chrN_start_end format). Junctions are pre-filtered to
those carried by at least one isoform with IF_overall > if_threshold in IsoformSwitchAnalyzeR.
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
    parser.add_argument("iso_file", help="IsoformSwitchAnalyzeR isoformFeatures.csv")
    parser.add_argument("-t", "--if_threshold", type=float, default=0.1,
                        help="Minimum IF_overall to consider an isoform expressed (default: 0.1)")
    parser.add_argument("-o", "--output", default="novel_coding_junction_sqtl_coloc_matches.tsv",
                        help="Output TSV file (default: novel_coding_junction_sqtl_coloc_matches.tsv)")
    args = parser.parse_args()

    iso = pl.read_csv(args.iso_file, null_values="NA")
    expressed_isoforms = set(
        iso.filter(pl.col("IF_overall") > args.if_threshold)["isoform_id"].unique().to_list()
    )
    print(f"Isoforms with IF_overall > {args.if_threshold}: {len(expressed_isoforms):,} of {iso['isoform_id'].n_unique():,}")

    sqtl = pl.read_csv(args.sqtl_matches, separator="\t", null_values="NA")

    def any_expressed(tx_ids: str) -> bool:
        return any(tx.strip() in expressed_isoforms for tx in tx_ids.split(";"))

    sqtl_filtered = sqtl.filter(
        pl.col("transcript_ids").map_elements(any_expressed, return_dtype=pl.Boolean)
    )
    print(f"sQTL junctions before expression filter: {sqtl['feature'].n_unique():,}")
    print(f"sQTL junctions after  expression filter: {sqtl_filtered['feature'].n_unique():,}")

    coloc = pl.read_csv(args.coloc_file, separator="\t", null_values="NA")
    result = sqtl_filtered.join(coloc, on="feature", how="inner")

    result.write_csv(args.output, separator="\t")
    print(f"Matched rows:      {len(result)}")
    print(f"Unique features:   {result['feature'].n_unique()}")
    print(f"Diseases present:  {sorted(result['disease'].unique().to_list())}")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
