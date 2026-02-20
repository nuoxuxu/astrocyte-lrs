#!/usr/bin/env python3
import re
import argparse

def get_transcript_ids(gtf_file):
    ids = set()
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m:
                ids.add(m.group(1))
    return ids

def get_lines_for_transcripts(gtf_file, transcript_ids):
    lines = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m and m.group(1) in transcript_ids:
                lines.append(line)
    return lines

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add non-coding transcripts from an original GTF to an ORFanage GTF. "
            "ORFanage only outputs transcripts predicted to be coding; this script "
            "identifies transcripts absent from the ORFanage output and appends them "
            "to produce a complete combined GTF."
        )
    )
    parser.add_argument(
        "-i", "--original",
        required=True,
        metavar="original.gtf",
        help="Original GTF file containing all transcripts (exons only, no CDS)."
    )
    parser.add_argument(
        "-r", "--orfanage",
        required=True,
        metavar="orfanage.gtf",
        help="ORFanage output GTF containing only coding transcripts with predicted CDS."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        metavar="combined.gtf",
        help="Output GTF file path for the combined result."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print summary statistics to stdout."
    )
    args = parser.parse_args()

    coding_ids    = get_transcript_ids(args.orfanage)
    all_ids       = get_transcript_ids(args.original)
    noncoding_ids = all_ids - coding_ids

    if args.verbose:
        print(f"Total transcripts in original GTF: {len(all_ids)}")
        print(f"Coding transcripts (ORFanage):     {len(coding_ids)}")
        print(f"Non-coding transcripts to add:     {len(noncoding_ids)}")

    noncoding_lines = get_lines_for_transcripts(args.original, noncoding_ids)

    with open(args.orfanage) as f:
        orfanage_lines = f.readlines()

    with open(args.output, "w") as out:
        out.writelines(orfanage_lines)
        out.writelines(noncoding_lines)

    if args.verbose:
        print(f"Combined GTF written to: {args.output}")

if __name__ == "__main__":
    main()