#!/usr/bin/env python
"""Side-by-side pie charts comparing cis-sQTL GWAS colocalization by trait:
   left = leafcutter_sqtl_coloc_matches.tsv, right = full BigBrain COLOC file."""

import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

LEAFCUTTER_FILE = "results/aim_2/leafcutter_sqtl_coloc_matches.tsv"
BIGBRAIN_FILE   = "data/BigBrain_cis_sQTL_COLOC_results_H4_0.5_with_LD.tsv"
OUT_FILE        = "figures/coloc_pie_by_trait.pdf"

COLORS = {
    "Height": "#4878D0",
    "SCZ":    "#EE854A",
    "RA":     "#6ACC65",
    "AD":     "#D65F5F",
    "MDD":    "#956CB4",
    "ALS":    "#8C613C",
    "PD":     "#DC7EC0",
    "BPD":    "#797979",
    "MS":     "#D5BB67",
}

ALL_TRAITS = ["Height", "SCZ", "RA", "AD", "MS", "ALS", "PD", "BPD", "MDD"]


def get_counts(path):
    df = pl.read_csv(path, separator="\t", null_values="NA")
    counts = {row["disease"]: row["n"] for row in df.group_by("disease").len().rename({"len": "n"}).to_dicts()}
    return counts


def make_pie(ax, counts, title):
    traits = [t for t in ALL_TRAITS if t in counts]
    ns     = [counts[t] for t in traits]
    total  = sum(ns)
    colors = [COLORS[t] for t in traits]

    ax.pie(
        ns,
        colors=colors,
        startangle=90,
        wedgeprops={"linewidth": 0.8, "edgecolor": "white"},
    )

    legend_labels = [f"{t}  (n={n}, {n/total*100:.1f}%)" for t, n in zip(traits, ns)]
    patches = [mpatches.Patch(color=c, label=l) for c, l in zip(colors, legend_labels)]
    ax.legend(
        handles=patches,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.30),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    ax.set_title(f"{title}\n(total = {total:,})", fontsize=11, pad=12)


lc_counts = get_counts(LEAFCUTTER_FILE)
bb_counts = get_counts(BIGBRAIN_FILE)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

make_pie(ax1, lc_counts, "Reactive astrocytes sQTLs")
make_pie(ax2, bb_counts, "All BigBrain sQTLs")

fig.suptitle("sQTL–GWAS colocalization by trait (PP.H4 ≥ 0.5)", fontsize=13, y=1.02)

plt.tight_layout()
plt.savefig(OUT_FILE, bbox_inches="tight")
print(f"Saved to {OUT_FILE}")
