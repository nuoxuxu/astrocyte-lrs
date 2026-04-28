#!/usr/bin/env python3
"""Filter formatted CPAT results to transcripts in the RiboTIE-filtered set.

CPAT is run on pre-RiboTIE transcripts (seq_ID = 'PB.10001.5').
RiboTIE-filtered FASTA has ORF-ID-style headers ('PB.10001.5_148').
This script:
  1. Builds a base_id -> [ribotie_ids] map from the filtered FASTA.
  2. Keeps only CPAT rows whose Sequence Name (base) has a RiboTIE ORF.
  3. For each matching RiboTIE ORF (a base may have >1), emits a row with
     Sequence Name replaced by the RiboTIE ID so IsoformSwitchAnalyzeR
     sees identifiers that match the switch list.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict


def parse_fasta_ids(fasta_path):
    """Return {base_id: [ribotie_id, ...]} mapping from FASTA headers."""
    mapping = defaultdict(list)
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            ribotie_id = line[1:].split()[0]
            base_id = re.sub(r"_\d+$", "", ribotie_id)
            mapping[base_id].append(ribotie_id)
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Filter CPAT results to RiboTIE-filtered transcripts"
    )
    parser.add_argument("cpat_results", help="Formatted CPAT file (convert_cpat_format.py output)")
    parser.add_argument("ribotie_fasta", help="RiboTIE-filtered FASTA (for ORF ID mapping)")
    parser.add_argument("-o", "--output", default="-", help="output file (default: stdout)")
    args = parser.parse_args()

    base_to_ribotie = parse_fasta_ids(args.ribotie_fasta)
    n_ribotie = sum(len(v) for v in base_to_ribotie.values())

    out = open(args.output, "w", newline="") if args.output != "-" else sys.stdout
    try:
        with open(args.cpat_results, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["Data ID", "Sequence Name", "RNA size", "ORF size",
                             "Ficket Score", "Hexamer Score", "Coding Probability", "Coding Label"])
            new_id = 0
            for row in reader:
                ribotie_ids = base_to_ribotie.get(row["Sequence Name"])
                if not ribotie_ids:
                    continue
                for ribotie_id in ribotie_ids:
                    writer.writerow([
                        new_id,
                        ribotie_id,
                        row["RNA size"],
                        row["ORF size"],
                        row["Ficket Score"],
                        row["Hexamer Score"],
                        row["Coding Probability"],
                        row["Coding Label"],
                    ])
                    new_id += 1
    finally:
        if args.output != "-":
            out.close()

    print(f"Wrote {new_id} rows (of {n_ribotie} RiboTIE transcripts)", file=sys.stderr)


if __name__ == "__main__":
    main()
