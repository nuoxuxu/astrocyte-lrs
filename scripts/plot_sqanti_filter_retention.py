#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

initial = 557807
counts = {
    "High\nStringency": 80589,
    "Mid\nStringency": 97337,
    "Low\nStringency": 162882,
}

labels = list(counts.keys())
kept = list(counts.values())
pcts = [n / initial * 100 for n in kept]

bar_height = 1.0
colors_kept = ["#2c6fad", "#2c6fad", "#2c6fad"]
colors_lost = ["#d0dce8", "#d0dce8", "#d0dce8"]

fig, ax = plt.subplots(figsize=(5, 4))

x = np.arange(len(labels))
width = 0.55

for i, (pct, n) in enumerate(zip(pcts, kept)):
    # Full bar (light = filtered out)
    ax.bar(x[i], bar_height, width=width, color=colors_lost[i], zorder=2)
    # Kept portion (dark, starts from bottom)
    ax.bar(x[i], pct / 100, width=width, color=colors_kept[i], zorder=3)
    # Label: percentage + count above bar
    ax.text(
        x[i],
        bar_height + 0.03,
        f"{pct:.1f}%\n(n={n:,})",
        ha="center",
        va="bottom",
        fontsize=9,
        fontweight="bold",
        color="#1a1a1a",
    )

ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylim(0, 1.35)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)
ax.set_ylabel("Proportion of initial transcripts", fontsize=10)
ax.set_title(
    f"Expression filter retention\n(initial set: {initial:,} transcripts)",
    fontsize=11,
    pad=10,
)

kept_patch = mpatches.Patch(color=colors_kept[0], label="Passed filter")
lost_patch = mpatches.Patch(color=colors_lost[0], label="Filtered out")
ax.legend(handles=[kept_patch, lost_patch], fontsize=9, loc="upper right", framealpha=0.9)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
out = "/scratch/nxu/astrocytes/figures/sqanti_filter_retention.pdf"
plt.savefig(out, bbox_inches="tight")
print(f"Saved to {out}")
