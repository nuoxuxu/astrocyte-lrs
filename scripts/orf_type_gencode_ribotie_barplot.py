#!/usr/bin/env python3
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

counts = Counter()
with open("nextflow_results/ribotie/gencode/merged/ribotie_res_merged.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        t = row.get("ORF_type", "").strip()
        if t:
            counts[t] += 1

order = [k for k, _ in counts.most_common()]
values = [counts[k] for k in order]

palette = [
    "#1565C0", "#42A5F5", "#90CAF9",
    "#FF6F00", "#FFA726",
    "#2E7D32", "#66BB6A",
    "#880E4F", "#CE93D8",
    "#795548", "#9E9E9E", "#BDBDBD",
]
colors = palette[:len(order)]

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.barh(order[::-1], values[::-1], color=colors[::-1], edgecolor="white", linewidth=0.5)

for bar, val in zip(bars, values[::-1]):
    ax.text(bar.get_width() + max(values) * 0.01, bar.get_y() + bar.get_height() / 2,
            f"{val:,}", va="center", ha="left", fontsize=9)

ax.set_xlabel("Number of ORFs", fontsize=11)
ax.set_title("ORF count by ORF_type (GENCODE RiboTIE)", fontsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(0, max(values) * 1.18)

plt.tight_layout()
out = "nextflow_results/quality/orf_type_gencode_ribotie_barplot.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
