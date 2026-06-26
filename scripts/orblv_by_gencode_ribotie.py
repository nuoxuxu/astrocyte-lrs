#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

df = pd.read_csv(
    "nextflow_results/quality/minlen_quality_metrics.tsv",
    sep="\t",
    usecols=["ORBLv", "in_gencode_ribotie"],
)
df = df[df["in_gencode_ribotie"].isin([True, False, "True", "False"])].copy()
df["ORBLv"] = pd.to_numeric(df["ORBLv"], errors="coerce")
df["in_gencode_ribotie"] = df["in_gencode_ribotie"].astype(str)
df = df.dropna(subset=["ORBLv"])

groups = {
    "True":  df[df["in_gencode_ribotie"] == "True"]["ORBLv"].values,
    "False": df[df["in_gencode_ribotie"] == "False"]["ORBLv"].values,
}

colors = {"True": "#2196F3", "False": "#E53935"}
labels = {
    "True":  f"In GENCODE RiboTIE (n={len(groups['True']):,})",
    "False": f"Not in GENCODE RiboTIE (n={len(groups['False']):,})",
}

fig, ax = plt.subplots(figsize=(7, 4.5))

x = np.linspace(0, 1, 500)
for key in ("True", "False"):
    vals = groups[key]
    # KDE via numpy histogram + smoothing
    counts, edges = np.histogram(vals, bins=100, range=(0, 1), density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    ax.fill_between(centers, counts, alpha=0.4, color=colors[key], linewidth=0)
    # Smooth KDE line using a simple Gaussian kernel via convolution
    sigma_bins = 3
    kernel = np.exp(-0.5 * (np.arange(-3*sigma_bins, 3*sigma_bins+1) / sigma_bins) ** 2)
    kernel /= kernel.sum()
    smoothed = np.convolve(counts, kernel, mode="same")
    ax.plot(centers, smoothed, color=colors[key], linewidth=1.8, label=labels[key])
    med = np.median(vals)
    ax.axvline(med, color=colors[key], linestyle="--", linewidth=1.0, alpha=0.8)

ax.set_xlabel("ORBLv score", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title("ORBLv score distribution by GENCODE RiboTIE membership", fontsize=12)
ax.set_xlim(0, 1)
ax.legend(fontsize=10, frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
out = "nextflow_results/quality/orblv_by_gencode_ribotie.png"
plt.savefig(out, dpi=150)
print(f"Saved {out}")

for key in ("True", "False"):
    vals = groups[key]
    print(f"{labels[key]}: median={np.median(vals):.4f}  mean={np.mean(vals):.4f}")
