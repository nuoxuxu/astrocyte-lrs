#!/usr/bin/env python3
"""
Convert hg38 chrom.sizes to use RefSeq-style accession names for non-primary contigs.

For each entry in chrom.sizes:
- Primary chromosomes (chr1-chr22, chrX, chrY, chrM): keep original name
- All others: replace name with the alias from chromAlias.txt that matches
  the pattern <alphanumeric>.number (no underscores)
"""

import re
import sys
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chromAlias", action='store', dest="chromAlias")
parser.add_argument("--chrom_sizes", action='store', dest="chrom_sizes")
args = parser.parse_args()

ALIAS_FILE = args.chromAlias
SIZES_FILE = args.chrom_sizes

PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}

# Matches e.g. "GL383545.1", "CM000663.2" — no underscores, ends with .digits
ACCESSION = re.compile(r"^[A-Za-z0-9]+\.[0-9]+$")


def load_alias_map(path):
    alias_map = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chr_name = parts[0]
            candidates = [x for x in parts[1:] if ACCESSION.match(x)]
            if candidates:
                alias_map[chr_name] = candidates[0]
    return alias_map


def main():
    alias_map = load_alias_map(ALIAS_FILE)

    with open(SIZES_FILE) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            chr_name, size = parts[0], parts[1]

            if chr_name in PRIMARY:
                out_name = chr_name
            else:
                out_name = alias_map.get(chr_name, chr_name)

            print(f"{out_name}\t{size}")


if __name__ == "__main__":
    main()
