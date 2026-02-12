#!/usr/bin/env python3
"""Subset GENCODE GTF to transcripts with class_code '=' in a gffcompare tmap file.

Extracts ref_id (GENCODE transcript IDs) from tmap rows where class_code == '=',
then filters the GENCODE GTF to keep only lines belonging to those transcripts
(transcript, exon, CDS, UTR, etc.) plus their parent gene lines.
"""

import argparse
import csv
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tmap", required=True, help="gffcompare .tmap file")
    parser.add_argument("--gtf", required=True, help="GENCODE annotation GTF")
    parser.add_argument("--output", required=True, help="Output GTF path")
    return parser.parse_args()


def get_matching_transcript_ids(tmap_path):
    """Return set of ref_id values where class_code == '='."""
    ids = set()
    with open(tmap_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["class_code"] == "=":
                ids.add(row["ref_id"])
    return ids


def extract_attribute(line, attr_name):
    """Extract an attribute value from a GTF attributes field."""
    match = re.search(rf'{attr_name} "([^"]+)"', line)
    return match.group(1) if match else None


def main():
    args = parse_args()

    transcript_ids = get_matching_transcript_ids(args.tmap)
    print(f"Found {len(transcript_ids)} transcript IDs with class_code '='", file=sys.stderr)

    # First pass: collect gene_ids that contain at least one matching transcript
    gene_ids = set()
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            tid = extract_attribute(line, "transcript_id")
            if tid in transcript_ids:
                gid = extract_attribute(line, "gene_id")
                if gid:
                    gene_ids.add(gid)

    # Second pass: write matching lines
    written = 0
    with open(args.gtf) as f, open(args.output, "w") as out:
        for line in f:
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.strip().split("\t")
            feature_type = fields[2]

            if feature_type == "gene":
                gid = extract_attribute(line, "gene_id")
                if gid in gene_ids:
                    out.write(line)
                    written += 1
            else:
                tid = extract_attribute(line, "transcript_id")
                if tid in transcript_ids:
                    out.write(line)
                    written += 1

    print(f"Wrote {written} GTF lines to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
