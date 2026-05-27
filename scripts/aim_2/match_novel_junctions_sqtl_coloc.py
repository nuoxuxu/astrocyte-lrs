#!/usr/bin/env python
"""
Cross-reference novel coding junction sQTLs with BigBrain pre-computed COLOC results.

Joins novel_coding_junction_sqtl_matches.tsv with BigBrain COLOC results (PP.H4 >= 0.5)
on the shared feature column (chrN_start_end format). Junctions are pre-filtered to
those carried by at least one isoform with IF_overall > 0.1 in IsoformSwitchAnalyzeR.
Output has one row per junction x disease x GWAS locus combination.
"""

import polars as pl

SQTL_FILE    = "results/aim_2/novel_coding_junction_sqtl_matches.tsv"
COLOC_FILE   = "data/BigBrain_cis_sQTL_COLOC_results_H4_0.5_with_LD.tsv"
ISO_FILE     = "nextflow_results/IsoformSwitchAnalyzeR/isoformFeatures.csv"
OUT_FILE     = "results/aim_2/novel_coding_junction_sqtl_coloc_matches.tsv"
IF_THRESHOLD = 0.1

# Build set of moderately expressed isoforms (IF_overall > threshold)
iso = pl.read_csv(ISO_FILE, null_values="NA")
expressed_isoforms = set(
    iso.filter(pl.col("IF_overall") > IF_THRESHOLD)["isoform_id"].unique().to_list()
)
print(f"Isoforms with IF_overall > {IF_THRESHOLD}: {len(expressed_isoforms):,} of {iso['isoform_id'].n_unique():,}")

# Keep junctions where at least one transcript_id passes the expression filter
sqtl = pl.read_csv(SQTL_FILE, separator="\t", null_values="NA")

def any_expressed(tx_ids: str) -> bool:
    return any(tx.strip() in expressed_isoforms for tx in tx_ids.split(";"))

sqtl_filtered = sqtl.filter(
    pl.col("transcript_ids").map_elements(any_expressed, return_dtype=pl.Boolean)
)
print(f"sQTL junctions before expression filter: {sqtl['feature'].n_unique():,}")
print(f"sQTL junctions after  expression filter: {sqtl_filtered['feature'].n_unique():,}")

coloc  = pl.read_csv(COLOC_FILE, separator="\t", null_values="NA")
result = sqtl_filtered.join(coloc, on="feature", how="inner")

result.write_csv(OUT_FILE, separator="\t")
print(f"Matched rows:      {len(result)}")
print(f"Unique features:   {result['feature'].n_unique()}")
print(f"Diseases present:  {sorted(result['disease'].unique().to_list())}")
print(f"Output written to: {OUT_FILE}")
