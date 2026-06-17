#!/usr/bin/env python3
"""
Fix malformed gene_name attributes in the collaborator GTF.

The bug: empty gene names are encoded as  gene_name ";
         instead of the correct           gene_name "";
         (the closing quote before the semicolon is missing).
"""

import argparse
import re
import sys

PATTERN = re.compile(r'gene_name ";')
FIX = 'gene_name "";'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_gtf")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    fixed = 0
    with open(args.input_gtf) as fin, open(args.output, "w") as fout:
        for line in fin:
            if PATTERN.search(line):
                line = PATTERN.sub(FIX, line)
                fixed += 1
            fout.write(line)

    print(f"Fixed {fixed:,} lines → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
