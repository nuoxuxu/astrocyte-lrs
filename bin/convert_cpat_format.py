#!/usr/bin/env python3
"""Convert CPAT.ORF_prob.best.tsv to IsoformSwitchAnalyzeR cpat_results.txt format."""

import argparse
import csv
import sys

HUMAN_CPAT_CUTOFF = 0.364


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="CPAT.ORF_prob.best.tsv")
    parser.add_argument("-o", "--output", default="-", help="output file (default: stdout)")
    parser.add_argument("--cutoff", type=float, default=HUMAN_CPAT_CUTOFF,
                        help=f"coding probability cutoff (default: {HUMAN_CPAT_CUTOFF})")
    args = parser.parse_args()

    out = open(args.output, "w", newline="") if args.output != "-" else sys.stdout

    try:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["Data ID", "Sequence Name", "RNA size", "ORF size",
                         "Ficket Score", "Hexamer Score", "Coding Probability", "Coding Label"])

        with open(args.input, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for i, row in enumerate(reader):
                prob = float(row["Coding_prob"])
                writer.writerow([
                    i,
                    row["seq_ID"],
                    row["mRNA"],
                    row["ORF"],
                    row["Fickett"],
                    row["Hexamer"],
                    prob,
                    "yes" if prob >= args.cutoff else "no",
                ])
    finally:
        if args.output != "-":
            out.close()


if __name__ == "__main__":
    main()
