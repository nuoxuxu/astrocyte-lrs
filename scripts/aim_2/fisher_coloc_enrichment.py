#!/usr/bin/env python
"""Fisher's exact test: are leafcutter sQTLs enriched for GWAS colocalization per trait?

2x2 table for each trait:
                      Colocalized   Not colocalized
  Leafcutter sQTLs        a              b
  Other BigBrain sQTLs    c              d

Universe = all unique features in BigBrain_cis_sQTL_ALL_top_assoc.tsv (N = 70,701)
Leafcutter set  = unique features in leafcutter_sqtl_matches.tsv (n = 1,239)
"""

import math
import polars as pl


def hypergeom_sf(k, N, K, n):
    """P(X >= k) where X ~ Hypergeometric(N, K, n). One-tailed p-value (greater)."""
    total = math.comb(N, n)
    return sum(
        math.comb(K, x) * math.comb(N - K, n - x)
        for x in range(k, min(K, n) + 1)
    ) / total


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    fdr = [0.0] * m
    for rank, i in enumerate(order, 1):
        fdr[i] = min(pvals[i] * m / rank, 1.0)
    # enforce monotonicity from right
    min_fdr = 1.0
    for i in reversed(order):
        min_fdr = min(fdr[i], min_fdr)
        fdr[i] = min_fdr
    return fdr

TOP_FILE      = "data/BigBrain_cis_sQTL_ALL_top_assoc.tsv"
LC_FILE       = "results/aim_2/leafcutter_sqtl_matches.tsv"
COLOC_FILE    = "data/BigBrain_cis_sQTL_COLOC_results_H4_0.5_with_LD.tsv"
LC_COLOC_FILE = "results/aim_2/leafcutter_sqtl_coloc_matches.tsv"

universe   = pl.read_csv(TOP_FILE,      separator="\t", null_values="NA")["feature"].unique()
lc_feats   = pl.read_csv(LC_FILE,       separator="\t", null_values="NA")["feature"].unique()
coloc      = pl.read_csv(COLOC_FILE,    separator="\t", null_values="NA")
lc_coloc   = pl.read_csv(LC_COLOC_FILE, separator="\t", null_values="NA")

N    = universe.len()   # 70,701
n_lc = lc_feats.len()  # 1,239

ALL_TRAITS = sorted(coloc["disease"].unique().to_list())

rows = []
for trait in ALL_TRAITS:
    a = lc_coloc.filter(pl.col("disease") == trait)["feature"].n_unique()
    b = n_lc - a

    c_total = coloc.filter(pl.col("disease") == trait)["feature"].n_unique()
    c = c_total - a                  # colocalized, not in leafcutter
    d = (N - n_lc) - c              # not colocalized, not in leafcutter

    # odds ratio
    odds_ratio = (a * d) / (b * c) if b * c > 0 else float("inf")
    # one-tailed p-value: P(X >= a) under hypergeometric null
    # N=total, K=total colocalized, n=leafcutter set size
    p_value = hypergeom_sf(a, N, a + c, n_lc)
    rows.append({"trait": trait, "a": a, "b": b, "c": c, "d": d,
                 "OR": odds_ratio, "p": p_value})

result = pl.DataFrame(rows)

fdr = bh_fdr(result["p"].to_list())
result = result.with_columns(pl.Series("FDR", fdr))

result = result.sort("p")

print(f"Universe (N):        {N:,}")
print(f"Leafcutter sQTLs:    {n_lc:,}\n")
print(result.select(["trait", "a", "b", "c", "d", "OR", "p", "FDR"]))
