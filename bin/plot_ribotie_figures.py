#!/usr/bin/env python3
"""Plot RiboTIE prediction proportions and expression density figures.

Produces two 3-panel figures:
  1. Proportion of transcripts with RiboTIE predictions by average expression level
  2. Distribution of average read counts by RiboTIE detection status
"""

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BINS = [0, 1, 5, 10, 50, 100, np.inf]
BIN_LABELS = ["(0,1]", "(1,5]", "(5,10]", "(10,50]", "(50,100]", ">100"]


# --------------- data loading helpers ---------------

def load_ribotie_union(stim_csv, unstim_csv):
    stim = pd.read_csv(stim_csv)
    unstim = pd.read_csv(unstim_csv)
    return set(stim["transcript_id"].unique()) | set(unstim["transcript_id"].unique())


def load_gtf_transcript_ids(gtf_path):
    tx_ids = set()
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] == "transcript":
                m = re.search(r'transcript_id "([^"]+)"', fields[8])
                if m:
                    tx_ids.add(m.group(1))
    return tx_ids


def avg_expression(expression_df):
    sample_cols = [c for c in expression_df.columns if c != "tname"]
    expression_df = expression_df.copy()
    expression_df["avg_count"] = expression_df[sample_cols].mean(axis=1)
    return expression_df


# --------------- proportion figure helpers ---------------

def compute_proportion(ribotie_tx, expression_df):
    expression_df = avg_expression(expression_df)
    df = expression_df[["tname", "avg_count"]].copy()
    df["has_ribotie"] = df["tname"].isin(ribotie_tx)
    df["count_bin"] = pd.cut(df["avg_count"], bins=BINS, labels=BIN_LABELS, right=True)

    grouped = df.groupby("count_bin", observed=False).agg(
        total=("has_ribotie", "size"),
        with_ribotie=("has_ribotie", "sum"),
    )
    grouped["proportion"] = grouped["with_ribotie"] / grouped["total"]
    grouped = grouped[grouped["total"] > 0]
    return grouped


def plot_proportion_panel(ax, grouped, title):
    tick_labels = [str(idx) for idx in grouped.index]
    ax.bar(range(len(grouped)), grouped["proportion"], color="#4C72B0", edgecolor="white")
    for i, (idx, row) in enumerate(grouped.iterrows()):
        ax.text(i, row["proportion"] + 0.005, f"n={int(row['total'])}",
                ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(grouped)))
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Average read count across samples")
    ax.set_ylabel("Proportion with RiboTIE prediction")
    ax.set_title(title)
    ax.set_ylim(0, min(grouped["proportion"].max() * 1.15, 1.0))


# --------------- density figure helpers ---------------

def prepare_density_df(ribotie_tx, expression_df):
    expression_df = avg_expression(expression_df)
    df = expression_df[["tname", "avg_count"]].copy()
    df["has_ribotie"] = df["tname"].isin(ribotie_tx)
    return df


def plot_density_panel(ax, df, title):
    detected = df.loc[df["has_ribotie"], "avg_count"]
    undetected = df.loc[~df["has_ribotie"], "avg_count"]

    detected_log = np.log10(detected[detected > 0])
    undetected_log = np.log10(undetected[undetected > 0])

    hist_bins = np.linspace(-1, 4, 60)

    ax.hist(undetected_log, bins=hist_bins, alpha=0.5, color="#D65F5F",
            label=f"No RiboTIE (n={len(undetected):,})")
    ax.hist(detected_log, bins=hist_bins, alpha=0.5, color="#4C72B0",
            label=f"RiboTIE detected (n={len(detected):,})")

    ax.set_xlabel("Average read count (log10)")
    ax.set_ylabel("Number of transcripts")
    ax.set_title(title)
    ax.legend(fontsize=9)


# --------------- per-dataset loaders ---------------

def load_isoseq_dataset(name, stim_csv, unstim_csv, orfanage_gtf, expression_parquet):
    ribotie_tx = load_ribotie_union(stim_csv, unstim_csv)
    orfanage_tx = load_gtf_transcript_ids(orfanage_gtf)
    expr = pd.read_parquet(expression_parquet)
    expr = expr[expr["tname"].isin(orfanage_tx)]
    grouped = compute_proportion(ribotie_tx, expr)
    density_df = prepare_density_df(ribotie_tx, expr)
    return {"name": name, "grouped": grouped, "density_df": density_df}


def load_gencode_dataset(name, stim_csv, unstim_csv, classification_parquet, expression_parquet):
    ribotie_tx = load_ribotie_union(stim_csv, unstim_csv)

    clf = pd.read_parquet(classification_parquet)
    isoform_to_gencode = clf[["isoform", "associated_transcript"]].copy()
    isoform_to_gencode = isoform_to_gencode[
        isoform_to_gencode["associated_transcript"].notna()
        & (isoform_to_gencode["associated_transcript"] != "novel")
    ]

    expr = pd.read_parquet(expression_parquet)
    sample_cols = [c for c in expr.columns if c != "tname"]

    gencode_expr = expr.merge(
        isoform_to_gencode, left_on="tname", right_on="isoform", how="inner"
    )
    gencode_agg = gencode_expr.groupby("associated_transcript")[sample_cols].sum().reset_index()
    gencode_agg.rename(columns={"associated_transcript": "tname"}, inplace=True)

    grouped = compute_proportion(ribotie_tx, gencode_agg)
    density_df = prepare_density_df(ribotie_tx, gencode_agg)
    return {"name": name, "grouped": grouped, "density_df": density_df}


# --------------- main ---------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--isoseq", nargs=5, action="append", default=[],
        metavar=("NAME", "STIM_CSV", "UNSTIM_CSV", "ORFANAGE_GTF", "EXPRESSION"),
        help="Isoseq dataset: label stim_csv unstim_csv orfanage_gtf expression_parquet",
    )
    parser.add_argument(
        "--gencode", nargs=5,
        metavar=("NAME", "STIM_CSV", "UNSTIM_CSV", "CLASSIFICATION", "EXPRESSION"),
        help="GENCODE dataset: label stim_csv unstim_csv classification_parquet expression_parquet",
    )
    parser.add_argument("--output_dir", default=".", help="Output directory for figures")
    args = parser.parse_args()

    datasets = []
    for iso in args.isoseq:
        name, stim, unstim, orf, expr = iso
        datasets.append(load_isoseq_dataset(name, stim, unstim, orf, expr))
    if args.gencode:
        name, stim, unstim, clf, expr = args.gencode
        datasets.append(load_gencode_dataset(name, stim, unstim, clf, expr))

    n = len(datasets)
    os.makedirs(args.output_dir, exist_ok=True)

    # --- Proportion figure ---
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6))
    if n == 1:
        axes = [axes]
    for ax, ds in zip(axes, datasets):
        plot_proportion_panel(ax, ds["grouped"], ds["name"])
    fig.suptitle(
        "Proportion of transcripts with RiboTIE predictions\nby average expression level",
        fontsize=14, y=1.02,
    )
    plt.tight_layout()
    fig.savefig(
        os.path.join(args.output_dir, "ribotie_proportion_by_expression.png"),
        dpi=200, bbox_inches="tight",
    )
    plt.close(fig)

    # --- Density figure ---
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 6))
    if n == 1:
        axes = [axes]
    for ax, ds in zip(axes, datasets):
        plot_density_panel(ax, ds["density_df"], ds["name"])
    fig.suptitle(
        "Distribution of average read counts\nby RiboTIE detection status",
        fontsize=14, y=1.02,
    )
    plt.tight_layout()
    fig.savefig(
        os.path.join(args.output_dir, "ribotie_expression_density.png"),
        dpi=200, bbox_inches="tight",
    )
    plt.close(fig)

    # --- Print summaries ---
    for ds in datasets:
        print(f"=== {ds['name']} ===")
        print(ds["grouped"])
        print()


if __name__ == "__main__":
    main()
