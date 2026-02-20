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
            "Restore transcripts removed by AGAT (which deduplicates transcripts with "
            "identical CDS) back to the output GTF. Takes the AGAT input and output GTFs, "
            "identifies missing transcript IDs, and appends their features from the input."
        )
    )
    parser.add_argument(
        "-i", "--agat-input",
        required=True,
        metavar="complete_orfanage.gtf",
        help="GTF file passed as input to AGAT (before deduplication)."
    )
    parser.add_argument(
        "-r", "--agat-output",
        required=True,
        metavar="agat_output.gtf",
        help="GTF file produced by AGAT (after deduplication and format conversion)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        metavar="restored.gtf",
        help="Output GTF file with removed transcripts added back."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print summary statistics to stdout."
    )
    args = parser.parse_args()

    pre_agat_ids  = get_transcript_ids(args.agat_input)
    post_agat_ids = get_transcript_ids(args.agat_output)
    removed_ids   = pre_agat_ids - post_agat_ids

    if args.verbose:
        print(f"Transcripts before AGAT: {len(pre_agat_ids)}")
        print(f"Transcripts after AGAT:  {len(post_agat_ids)}")
        print(f"Transcripts to restore:  {len(removed_ids)}")

    restored_lines = get_lines_for_transcripts(args.agat_input, removed_ids)

    with open(args.agat_output) as f:
        agat_lines = f.readlines()

    with open(args.output, "w") as out:
        out.writelines(agat_lines)
        out.writelines(restored_lines)

    if args.verbose:
        print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
