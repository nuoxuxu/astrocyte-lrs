#!/usr/bin/env python
"""
Cross-reference significant leafcutter sQTLs with BigBrain pre-computed COLOC results.

Joins leafcutter_sqtl_matches.tsv with BigBrain COLOC results (PP.H4 >= 0.5)
on the shared feature column (chrN_start_end format). Output has one row per
intron cluster x disease x GWAS locus combination.
"""

import polars as pl

SQTL_FILE  = "results/aim_2/leafcutter_sqtl_matches.tsv"
COLOC_FILE = "data/BigBrain_cis_sQTL_COLOC_results_H4_0.5_with_LD.tsv"
OUT_FILE   = "results/aim_2/leafcutter_sqtl_coloc_matches.tsv"

sqtl  = pl.read_csv(SQTL_FILE,  separator="\t", null_values="NA")
coloc = pl.read_csv(COLOC_FILE, separator="\t", null_values="NA")

result = sqtl.join(coloc, on="feature", how="inner")

result.write_csv(OUT_FILE, separator="\t")
print(f"sQTL rows:         {len(sqtl)}")
print(f"COLOC rows:        {len(coloc)}")
print(f"Matched rows:      {len(result)}")
print(f"Unique features:   {result['feature'].n_unique()}")
print(f"Diseases present:  {sorted(result['disease'].unique().to_list())}")
print(f"Output written to: {OUT_FILE}")
