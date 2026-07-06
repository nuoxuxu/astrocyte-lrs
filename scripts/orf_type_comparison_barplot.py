#!/usr/bin/env python3
"""Grouped bar chart comparing ORF_type_ORFanage vs ORF_type_RiboTIE counts."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

orfanage_counts = Counter()
gencode_counts = Counter()

with open("nextflow_results/quality/no_minlen_orf_type_gencode.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        orfanage_counts[row["ORF_type_ORFanage"].strip()] += 1
        gencode_counts[row["ORF_type_RiboTIE"].strip()] += 1

# All categories union, ordered: annotated CDS first, then by ORFanage count desc,
# then no_template_* at the end
CANONICAL_ORDER = [
    "annotated CDS",
    "C-terminal extension",
    "C-terminal truncation",
    "N-terminal extension",
    "N-terminal truncation",
    "uORF",
    "uoORF",
    "dORF",
    "doORF",
    "intORF",
    "lncRNA-ORF",
    "varRNA-ORF",
]
all_cats = set(orfanage_counts) | set(gencode_counts)
no_template = sorted(c for c in all_cats if c.startswith("no_template_"))
extra = sorted(c for c in all_cats if c not in CANONICAL_ORDER and not c.startswith("no_template_"))
order = [c for c in CANONICAL_ORDER if c in all_cats] + extra + no_template

orfanage_vals = [orfanage_counts.get(c, 0) for c in order]
gencode_vals  = [gencode_counts.get(c, 0) for c in order]

# Horizontal grouped bars
n = len(order)
y = np.arange(n)
height = 0.36

fig, ax = plt.subplots(figsize=(9, 6.5))

bars_o = ax.barh(y + height / 2, orfanage_vals, height, label="ORF_type_ORFanage",
                 color="#1565C0", alpha=0.88)
bars_g = ax.barh(y - height / 2, gencode_vals,  height, label="ORF_type_RiboTIE",
                 color="#E65100", alpha=0.88)

# Value labels
for bar, val in zip(bars_o, orfanage_vals):
    if val > 0:
        ax.text(bar.get_width() + 80, bar.get_y() + bar.get_height() / 2,
                f"{val:,}", va="center", ha="left", fontsize=7.5)
for bar, val in zip(bars_g, gencode_vals):
    if val > 0:
        ax.text(bar.get_width() + 80, bar.get_y() + bar.get_height() / 2,
                f"{val:,}", va="center", ha="left", fontsize=7.5)

ax.set_yticks(y)
ax.set_yticklabels(order, fontsize=9)
ax.set_xlabel("Number of ORFs", fontsize=11)
ax.set_title("ORF count by type — ORFanage vs GENCODE (no_minlen)", fontsize=12)
ax.legend(fontsize=9, loc="lower right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
max_val = max(max(orfanage_vals), max(gencode_vals))
ax.set_xlim(0, max_val * 1.20)

plt.tight_layout()
out = "nextflow_results/quality/orf_type_comparison_barplot_no_minlen.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
