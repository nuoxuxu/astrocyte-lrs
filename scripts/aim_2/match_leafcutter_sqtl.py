#!/usr/bin/env python
"""
Match significant leafcutter intron clusters with BigBrain cis-sQTLs.

Steps:
  1. Filter leafcutter_ds_cluster_significance.txt: status=="Success" and p.adjust<0.05
  2. Extract cluster names (strip chr prefix), match against leafcutter.clu2gene.txt col 1
  3. Subtract 1 from end coordinate; build feature string chr_start_end
  4. Join with BigBrain_cis_sQTL_ALL_top_assoc.tsv on feature; keep meta-analysis columns
  5. Write to results/aim_2/leafcutter_sqtl_matches.tsv
"""

import polars as pl

SIG_FILE  = "/scratch/nxu/leafcutter/example_data/leafcutter_ds_cluster_significance.txt"
CLU2GENE  = "/scratch/nxu/leafcutter/example_data/tmpdir_astrocytes/leafcutter.clu2gene.txt"
SQTL_FILE = "data/BigBrain_cis_sQTL_ALL_top_assoc.tsv"
OUT_FILE  = "results/aim_2/leafcutter_sqtl_matches.tsv"

META_COLS = [
    "feature", "gene_ids", "gene_names",
    "variant_id", "chr", "pos", "ref", "alt", "Allele",
    "fixed_beta", "fixed_sd", "fixed_z", "Random_Z",
    "Fixed_P", "Random_P", "Fixed_bonf", "Random_bonf",
    "Fixed_FDR", "Random_FDR", "Bonf_FDR",
]

# Step 1: filter significant clusters
sig = pl.read_csv(SIG_FILE, separator="\t", null_values="NA")
sig_filtered = sig.filter(
    (pl.col("status") == "Success") & (pl.col("p.adjust") < 0.05)
)
cluster_names = (
    sig_filtered["cluster"]
    .str.split(":")
    .list.get(1)
    .to_list()
)
print(f"Significant clusters (p.adjust < 0.05): {len(cluster_names)}")

# Step 2: filter clu2gene to significant clusters
clu2gene = pl.read_csv(
    CLU2GENE,
    separator=" ",
    has_header=False,
    new_columns=["clu_name", "chr", "start", "end", "gene_id", "m1", "m2", "m3", "m4"],
)
matched = clu2gene.filter(pl.col("clu_name").is_in(cluster_names))
print(f"clu2gene rows matched: {len(matched)}")

# Step 3: build feature string (subtract 1 from end coordinate)
matched = matched.with_columns((pl.col("end") - 1).alias("end"))
matched = matched.with_columns(
    (
        pl.col("chr") + "_"
        + pl.col("start").cast(pl.String) + "_"
        + pl.col("end").cast(pl.String)
    ).alias("feature")
)
features = set(matched["feature"].to_list())

# Step 4: load sQTL file and filter to matching features
sqtl = pl.read_csv(SQTL_FILE, separator="\t", null_values="NA")
sqtl_filtered = sqtl.filter(pl.col("feature").is_in(features))
keep_cols = [c for c in META_COLS if c in sqtl_filtered.columns]
result = sqtl_filtered.select(keep_cols)
print(f"sQTL rows matched: {len(result)}")

# Step 5: write output
result.write_csv(OUT_FILE, separator="\t")
print(f"Output written to: {OUT_FILE}")
