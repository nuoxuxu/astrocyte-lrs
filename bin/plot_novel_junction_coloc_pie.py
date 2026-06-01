#!/usr/bin/env python3
"""Side-by-side pie charts comparing cis-sQTL GWAS colocalization by trait:
   left = novel coding junctions intersected with COLOC results,
   right = full BigBrain COLOC file."""

import argparse

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import polars as pl

COLORS = {
    "SCZ":    "#EE854A",
    "RA":     "#6ACC65",
    "AD":     "#D65F5F",
    "MDD":    "#956CB4",
    "ALS":    "#8C613C",
    "PD":     "#DC7EC0",
    "BPD":    "#797979",
    "MS":     "#D5BB67",
}

ALL_TRAITS = ["SCZ", "RA", "AD", "MS", "ALS", "PD", "BPD", "MDD"]


def get_counts(df_or_path):
    df = pl.read_csv(df_or_path, separator="\t", null_values="NA") if isinstance(df_or_path, str) else df_or_path
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


def main():
    parser = argparse.ArgumentParser(
        description="Side-by-side pie charts comparing sQTL GWAS colocalization by trait."
    )
    parser.add_argument("novel_coloc_file", help="Novel coding junction sQTL COLOC matches TSV")
    parser.add_argument("bigbrain_file", help="Full BigBrain COLOC results TSV (PP.H4 >= 0.5)")
    parser.add_argument("-o", "--output", default="novel_coding_junction_coloc_pie.pdf",
                        help="Output PDF file (default: novel_coding_junction_coloc_pie.pdf)")
    args = parser.parse_args()

    novel_coloc = pl.read_csv(args.novel_coloc_file, separator="\t", null_values="NA")
    print(f"Novel coding junction COLOC rows:   {len(novel_coloc):,}")

    novel_counts = get_counts(novel_coloc)
    bb_counts    = get_counts(args.bigbrain_file)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    make_pie(ax1, novel_counts, "Novel coding junction sQTLs")
    make_pie(ax2, bb_counts,    "All BigBrain sQTLs")

    fig.suptitle("sQTL–GWAS colocalization by trait (PP.H4 ≥ 0.5)", fontsize=13, y=1.02)

    plt.tight_layout()
    plt.savefig(args.output, bbox_inches="tight")
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main()
