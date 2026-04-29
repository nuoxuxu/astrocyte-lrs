#!/usr/bin/env python3
"""Plot overlaid CDS length distributions for RiboTIE ORF predictions across orfanage modes."""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from utils import read_fasta

COLORS = ["#4C72B0", "#D65F5F"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastas", nargs="+", required=True)
    parser.add_argument("--modes", nargs="+", required=True)
    parser.add_argument("--output_png", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    assert len(args.fastas) == len(args.modes), "--fastas and --modes must have the same length"

    fig, ax = plt.subplots(figsize=(8, 5))

    for fasta, mode, color in zip(args.fastas, args.modes, COLORS):
        df = read_fasta(fasta)
        cds_lengths = [len(seq) * 3 for seq in df["seq"].to_list()]
        ax.hist(cds_lengths, bins=50, alpha=0.6, color=color, edgecolor="none", label=mode)

    ax.set_xlabel("CDS Length (nt)")
    ax.set_ylabel("Count")
    ax.set_title("RiboTIE ORF CDS Length Distribution")
    ax.legend()
    plt.tight_layout()
    plt.savefig(args.output_png, dpi=200)
    plt.close()


if __name__ == "__main__":
    main()
