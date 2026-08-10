#!/usr/bin/env python3
"""
Fix malformed gene_name attributes and remove redundant duplicate ORFs from collaborator GTF.

The bug: empty gene names are encoded as  gene_name ";
         instead of the correct           gene_name "";
         (the closing quote before the semicolon is missing).

The duplication: filtered_output.gtf contains ~20-27 copies of each ORF feature with identical
coordinates but different ribotie_score values, likely from concatenating pre-deduplicated
RiboTIE predictions. This script keeps only the first occurrence of each unique feature.
"""

import argparse
import re
import sys

PATTERN = re.compile(r'gene_name ";')
FIX = 'gene_name "";'


def get_transcript_id(line):
    m = re.search(r'transcript_id "([^"]+)"', line)
    return m.group(1) if m else None


def deduplicate_gtf_lines(lines):
    """Deduplicate GTF lines by keeping first occurrence of each (chrom, feature, start, end, transcript_id)."""
    seen = set()
    deduped = []
    for line in lines:
        if line.startswith('#'):
            deduped.append(line)
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 8:
            deduped.append(line)
            continue
        chrom, source, feature, start, end, score, strand, frame = parts[:8]
        tid = get_transcript_id(line)
        key = (chrom, feature, start, end, tid)
        if key not in seen:
            seen.add(key)
            deduped.append(line)
    return deduped


def verify_no_duplicates(gtf_file):
    """Verify that the output GTF has no duplicate features per transcript."""
    seen = {}
    duplicates = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom, source, feature, start, end, score, strand, frame = parts[:8]
            tid = get_transcript_id(line)
            key = (chrom, feature, start, end, tid)
            if key in seen:
                duplicates.append((key, seen[key], line.rstrip()))
            else:
                seen[key] = line.rstrip()
    return duplicates


def main():
    parser = argparse.ArgumentParser(
        description="Fix gene_name formatting and remove redundant duplicate ORFs from GTF"
    )
    parser.add_argument("input_gtf")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    fixed_names = 0
    all_lines = []
    with open(args.input_gtf) as fin:
        for line in fin:
            if PATTERN.search(line):
                line = PATTERN.sub(FIX, line)
                fixed_names += 1
            all_lines.append(line)

    deduped_lines = deduplicate_gtf_lines(all_lines)
    removed_duplicates = len(all_lines) - len(deduped_lines)

    with open(args.output, "w") as fout:
        for line in deduped_lines:
            fout.write(line)

    print(f"Fixed {fixed_names:,} gene_name formatting issues", file=sys.stderr)
    print(f"Removed {removed_duplicates:,} duplicate feature lines", file=sys.stderr)
    print(f"Output → {args.output}", file=sys.stderr)

    # Verify deduplication worked
    duplicates = verify_no_duplicates(args.output)
    if duplicates:
        print(f"\nERROR: Found {len(duplicates)} duplicate features in output GTF:", file=sys.stderr)
        for key, first_line, second_line in duplicates[:5]:
            chrom, feature, start, end, tid = key
            print(f"  {chrom}:{feature}:{start}-{end} (transcript_id={tid})", file=sys.stderr)
        if len(duplicates) > 5:
            print(f"  ... and {len(duplicates) - 5} more", file=sys.stderr)
        sys.exit(1)

    print(f"✓ Verification passed: no duplicate features found", file=sys.stderr)


if __name__ == "__main__":
    main()
