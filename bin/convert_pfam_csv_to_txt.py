#!/usr/bin/env python3
"""Convert pfam_scan.py CSV output to the pfam_scan.pl whitespace-separated
format expected by IsoformSwitchAnalyzeR::analyzePFAM.

Usage:
    python convert_pfam_csv_to_txt.py input.csv output.txt
"""

import argparse
import csv
import sys


HEADER_COMMENT = """\
# pfam_scan.pl-compatible format converted from pfam_scan.py CSV output
#
# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
"""


def main():
    parser = argparse.ArgumentParser(
        description="Convert pfam_scan.py CSV to pfam_scan.pl TXT format"
    )
    parser.add_argument("input", help="Input CSV file from pfam_scan.py")
    parser.add_argument("output", help="Output TXT file for analyzePFAM")
    args = parser.parse_args()

    with open(args.input) as fin, open(args.output, "w") as fout:
        fout.write(HEADER_COMMENT)
        reader = csv.DictReader(fin)
        for row in reader:
            clan = row["clan"].replace("No_Clan", "No_clan")
            fields = [
                row["seq_id"],
                row["aln_start"],
                row["aln_end"],
                row["env_start"],
                row["env_end"],
                row["hmm_acc"],
                row["hmm_name"],
                row["type"],
                row["hmm_start"],
                row["hmm_end"],
                row["hmm_length"],
                row["score"],
                row["evalue"],
                row["significance"],
                clan,
            ]
            # Join fields with spaces to match pfam_scan.pl format.
            # Use explicit space separation to avoid merging when
            # field values exceed their expected width.
            line = " ".join(fields)
            fout.write(line + "\n")


if __name__ == "__main__":
    main()
