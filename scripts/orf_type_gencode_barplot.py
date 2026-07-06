#!/usr/bin/env python3
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

counts = Counter()
with open("nextflow_results/quality/no_minlen_quality_metrics.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        t = row["ORF_type_RiboTIE"].strip()
        if t:
            counts[t] += 1

# Ordered by count descending
order = [k for k, _ in counts.most_common()]
values = [counts[k] for k in order]

# Color palette: one color per category group
palette = [
    "#1565C0",  # annotated CDS
    "#42A5F5",  # C-terminal extension
    "#90CAF9",  # C-terminal truncation
    "#FF6F00",  # N-terminal extension
    "#FFA726",  # N-terminal truncation
    "#2E7D32",  # uORF
    "#66BB6A",  # uoORF
    "#880E4F",  # dORF
    "#CE93D8",  # doORF
    "#795548",  # intORF
    "#9E9E9E",  # lncRNA-ORF
    "#BDBDBD",  # varRNA-ORF
]
colors = palette[:len(order)]

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.barh(order[::-1], values[::-1], color=colors[::-1], edgecolor="white", linewidth=0.5)

for bar, val in zip(bars, values[::-1]):
    ax.text(bar.get_width() + 150, bar.get_y() + bar.get_height() / 2,
            f"{val:,}", va="center", ha="left", fontsize=9)

ax.set_xlabel("Number of ORFs", fontsize=11)
ax.set_title("ORF count by ORF_type_RiboTIE (no_minlen)", fontsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xlim(0, max(values) * 1.18)

plt.tight_layout()
out = "nextflow_results/quality/orf_type_gencode_barplot_no_minlen.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
