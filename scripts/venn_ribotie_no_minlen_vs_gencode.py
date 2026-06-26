#!/usr/bin/env python3
"""Venn diagram: no_minlen vs gencode RiboTIE ORFs, matched by protein_seq."""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def read_protein_seqs(path):
    seqs = set()
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq = row.get("protein_seq", "").strip()
            if seq:
                seqs.add(seq)
    return seqs

no_minlen = read_protein_seqs(
    "nextflow_results/ribotie/no_minlen/merged/ribotie_res_merged.csv"
)
gencode = read_protein_seqs(
    "nextflow_results/ribotie/gencode/merged/ribotie_res_merged.csv"
)

only_no_minlen = len(no_minlen - gencode)
only_gencode   = len(gencode - no_minlen)
shared         = len(no_minlen & gencode)

print(f"no_minlen only : {only_no_minlen:,}")
print(f"gencode only   : {only_gencode:,}")
print(f"shared         : {shared:,}")
print(f"no_minlen total: {len(no_minlen):,}")
print(f"gencode total  : {len(gencode):,}")

# --- draw Venn diagram manually (no matplotlib_venn dependency) ---
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.set_aspect("equal")
ax.axis("off")

cx1, cx2, cy, r = 0.38, 0.62, 0.5, 0.28

c1 = mpatches.Circle((cx1, cy), r, transform=ax.transAxes,
                      color="#2196F3", alpha=0.45, zorder=2)
c2 = mpatches.Circle((cx2, cy), r, transform=ax.transAxes,
                      color="#E53935", alpha=0.45, zorder=2)
ax.add_patch(c1)
ax.add_patch(c2)

# labels inside circles
ax.text(cx1 - r * 0.52, cy, f"{only_no_minlen:,}",
        ha="center", va="center", fontsize=14, fontweight="bold",
        transform=ax.transAxes, zorder=3)
ax.text((cx1 + cx2) / 2, cy, f"{shared:,}",
        ha="center", va="center", fontsize=14, fontweight="bold",
        transform=ax.transAxes, zorder=3)
ax.text(cx2 + r * 0.52, cy, f"{only_gencode:,}",
        ha="center", va="center", fontsize=14, fontweight="bold",
        transform=ax.transAxes, zorder=3)

# set labels
ax.text(cx1 - r * 0.3, cy + r + 0.06, "no_minlen\nIso-Seq",
        ha="center", va="bottom", fontsize=11, color="#1565C0",
        transform=ax.transAxes)
ax.text(cx2 + r * 0.3, cy + r + 0.06, "GENCODE\nreference",
        ha="center", va="bottom", fontsize=11, color="#B71C1C",
        transform=ax.transAxes)

ax.set_title("RiboTIE ORF overlap (matched by protein sequence)",
             fontsize=12, pad=12)

plt.tight_layout()
out = "nextflow_results/quality/venn_ribotie_no_minlen_vs_gencode.png"
plt.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved {out}")
