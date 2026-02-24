#!/usr/bin/env python3
"""Combine GENCODE and novel transcripts into a single GTF for RiboTIE input.

For each transcript in final_classification.parquet:
  - If associated_transcript is a valid GENCODE transcript ID → use the GENCODE version
  - Otherwise (novel, NA, null) → use the transcript from final_transcripts.gtf

The total number of unique transcripts in the output should equal the number in
final_transcripts.gtf. A warning is printed if they differ (e.g. multiple isoforms
map to the same GENCODE transcript).
"""
import argparse
import re
import sys

import polars as pl


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final_transcripts_gtf", required=True,
                        help="final_transcripts.gtf from filter_by_expression")
    parser.add_argument("--final_classification", required=True,
                        help="final_classification.parquet")
    parser.add_argument("--annotation_gtf", required=True,
                        help="GENCODE annotation GTF")
    parser.add_argument("--output_gtf", required=True,
                        help="Output combined GTF path")
    return parser.parse_args()


def get_transcript_ids(gtf_path):
    """Return the set of all transcript_id values in a GTF file."""
    ids = set()
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m:
                ids.add(m.group(1))
    return ids


def extract_records(gtf_path, transcript_ids, include_gene_entries=False):
    """Return GTF lines whose transcript_id is in transcript_ids.

    If include_gene_entries is True, gene-level lines whose gene_id is a parent
    of any selected transcript are also included.
    """
    gene_ids: set[str] = set()

    if include_gene_entries:
        with open(gtf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                m_tid = re.search(r'transcript_id "([^"]+)"', line)
                if m_tid and m_tid.group(1) in transcript_ids:
                    m_gid = re.search(r'gene_id "([^"]+)"', line)
                    if m_gid:
                        gene_ids.add(m_gid.group(1))

    lines: list[str] = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            feature = fields[2]
            if feature == "gene":
                if include_gene_entries:
                    m = re.search(r'gene_id "([^"]+)"', line)
                    if m and m.group(1) in gene_ids:
                        lines.append(line)
            else:
                m = re.search(r'transcript_id "([^"]+)"', line)
                if m and m.group(1) in transcript_ids:
                    lines.append(line)
    return lines


def main():
    args = parse_args()

    cls = pl.read_parquet(args.final_classification)

    # Isoforms with a valid GENCODE associated_transcript
    has_gencode = (
        pl.col("associated_transcript").is_not_null()
        & ~pl.col("associated_transcript").is_in(["novel", "NA", ""])
    )
    matched_cls = cls.filter(has_gencode)
    novel_cls = cls.filter(~has_gencode)

    gencode_tx_ids: set[str] = set(matched_cls["associated_transcript"].unique().to_list())
    novel_isoform_ids: set[str] = set(novel_cls["isoform"].to_list())

    # Verify each associated_transcript actually exists in the annotation GTF
    annotation_tx_ids = get_transcript_ids(args.annotation_gtf)
    valid_gencode_ids = gencode_tx_ids & annotation_tx_ids
    invalid_ids = gencode_tx_ids - valid_gencode_ids
    if invalid_ids:
        print(
            f"WARNING: {len(invalid_ids)} associated_transcript ID(s) not found in "
            "annotation GTF — treating corresponding isoforms as novel.",
            file=sys.stderr,
        )
        extra_novel = matched_cls.filter(
            pl.col("associated_transcript").is_in(list(invalid_ids))
        )
        novel_isoform_ids.update(extra_novel["isoform"].to_list())

    print(f"Transcripts from GENCODE:           {len(valid_gencode_ids)}", file=sys.stderr)
    print(f"Novel transcripts (final_transcripts.gtf): {len(novel_isoform_ids)}", file=sys.stderr)

    # Extract records from each source
    # GENCODE entries include gene-level lines; final_transcripts.gtf typically
    # does not contain gene-level lines (isoseq collapse output), so we skip them.
    gencode_lines = extract_records(args.annotation_gtf, valid_gencode_ids, include_gene_entries=True)
    novel_lines = extract_records(args.final_transcripts_gtf, novel_isoform_ids, include_gene_entries=False)

    with open(args.output_gtf, "w") as out:
        for line in gencode_lines:
            out.write(line)
        for line in novel_lines:
            out.write(line)

    total = len(valid_gencode_ids) + len(novel_isoform_ids)
    n_final = len(get_transcript_ids(args.final_transcripts_gtf))
    print(f"Total unique transcripts written:   {total}", file=sys.stderr)
    print(f"Transcripts in final_transcripts.gtf: {n_final}", file=sys.stderr)
    if total != n_final:
        print(
            f"WARNING: transcript counts differ (output {total} vs final_transcripts.gtf {n_final}). "
            "This is expected if multiple isoforms map to the same GENCODE transcript.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
