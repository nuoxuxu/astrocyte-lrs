#!/usr/bin/env python3
"""
Generate GTF files containing ORFs from two Ribo-seq studies, both based on GENCODE v45.

Study 1: 41587_2022_1369_MOESM2_ESM.xlsx, sheet "S2. PHASE I Ribo-seq ORFs"
  - Coordinates are 1-based (GTF format), semicolon-separated for multi-exon ORFs

Study 2: media-1.xlsx, sheet "Supplementary Table 2"
  - Coordinates are 0-based (BED format), comma-separated for multi-exon ORFs
  - Converted to 1-based by: start_gtf = start_bed + 1, end_gtf = end_bed

Stop codon convention: CDS coordinates in the source data include the stop codon.
This script trims the terminal stop codon from CDS coordinates using the same
trim_stop_codon function used for Ribo-TISH GTFs, so the output follows the
standard GTF/GENCODE convention (stop codon excluded from CDS).
"""

import sys
from pathlib import Path

import pandas as pd
import polars as pl

sys.path.insert(0, str(Path(__file__).parent.parent))
from bin.trim_ribotish_stop_codon import trim_stop_codon

DATA_DIR = Path("data")
OUT_DIR = Path("data")

GTF_COLUMNS = ["seqname", "source", "feature", "start", "end",
               "score", "strand", "frame", "attributes"]


def make_attributes(**kwargs) -> str:
    """Format GTF attribute string from keyword arguments."""
    parts = []
    for key, value in kwargs.items():
        if value is not None and str(value) != "nan":
            parts.append(f'{key} "{value}"')
    return "; ".join(parts)


def calc_frames(exon_lengths: list[int], strand: str) -> list[int]:
    """
    Calculate CDS frame for each exon block.

    GTF frame: number of bases at the start of this feature that complete
    a codon begun in the previous feature.

    For + strand: process left to right (as given).
    For - strand: the 5' end is the rightmost exon, so process right to left
                  but the frame is recorded for each exon in genomic order.
    """
    if strand == "+":
        frames = []
        cumulative = 0
        for length in exon_lengths:
            frames.append((3 - (cumulative % 3)) % 3)
            cumulative += length
        return frames
    else:  # minus strand: 5' -> 3' goes right to left
        frames_reversed = []
        cumulative = 0
        for length in reversed(exon_lengths):
            frames_reversed.append((3 - (cumulative % 3)) % 3)
            cumulative += length
        # frames_reversed[i] corresponds to exon at index (n-1-i)
        return list(reversed(frames_reversed))


def parse_coords(starts_str: str, ends_str: str, sep: str) -> tuple[list[int], list[int]]:
    """Parse semicolon- or comma-separated coordinate strings into integer lists."""
    starts = [int(s.strip()) for s in str(starts_str).split(sep)]
    ends = [int(e.strip()) for e in str(ends_str).split(sep)]
    assert len(starts) == len(ends), f"Mismatched starts/ends: {starts_str} | {ends_str}"
    return starts, ends


def build_orf_records(
    chrom: str,
    strand: str,
    starts_1based: list[int],
    ends_1based: list[int],
    orf_id: str,
    gene_id: str,
    gene_name: str,
    transcript_id: str,
    orf_type: str,
    gene_biotype: str,
    source: str,
) -> list[dict]:
    """Return a list of GTF row dicts (transcript + CDS lines) for one ORF."""
    tx_start = min(starts_1based)
    tx_end = max(ends_1based)

    base_attrs = make_attributes(
        gene_id=gene_id,
        gene_name=gene_name,
        transcript_id=transcript_id,
        orf_id=orf_id,
        orf_type=orf_type,
        gene_biotype=gene_biotype,
    )

    rows = []

    # Transcript-level record
    rows.append({
        "seqname": chrom, "source": source, "feature": "transcript",
        "start": tx_start, "end": tx_end, "score": ".", "strand": strand,
        "frame": ".", "attributes": base_attrs,
    })

    # CDS records (one per exon block)
    exon_lengths = [e - s + 1 for s, e in zip(starts_1based, ends_1based)]
    frames = calc_frames(exon_lengths, strand)

    for start, end, frame in zip(starts_1based, ends_1based, frames):
        rows.append({
            "seqname": chrom, "source": source, "feature": "CDS",
            "start": start, "end": end, "score": ".", "strand": strand,
            "frame": str(frame), "attributes": base_attrs,
        })

    return rows


# ---------------------------------------------------------------------------
# Study 1
# ---------------------------------------------------------------------------

def process_study1(xlsx_path: Path, out_path: Path) -> None:
    df = pd.read_excel(xlsx_path, sheet_name="S2. PHASE I Ribo-seq ORFs", header=0)
    print(f"Study 1: loaded {len(df)} ORFs")

    rows: list[dict] = []
    for _, row in df.iterrows():
        chrom = f"chr{row['chrm']}"
        strand = str(row["strand"])
        starts, ends = parse_coords(row["starts"], row["ends"], sep=";")
        rows.extend(build_orf_records(
            chrom=chrom, strand=strand,
            starts_1based=starts, ends_1based=ends,
            orf_id=str(row["orf_name"]),
            gene_id=str(row.get("gene_id", "")),
            gene_name=str(row.get("gene_name", "")),
            transcript_id=str(row.get("transcript", "")),
            orf_type=str(row.get("orf_biotype", "")),
            gene_biotype=str(row.get("gene_biotype", "")),
            source="study1_riboseq",
        ))

    gtf_df = pl.DataFrame(rows, schema={c: pl.String for c in GTF_COLUMNS})
    gtf_df = gtf_df.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.col("attributes").str.extract(r'orf_id "([^"]+)"').alias("orf_id"),
    )

    n_before = gtf_df.filter(pl.col("feature") == "CDS").height
    trimmed = trim_stop_codon(gtf_df, group_key="orf_id").drop("orf_id")
    n_after = trimmed.filter(pl.col("feature") == "CDS").height
    dropped = n_before - n_after
    print(f"Study 1: CDS rows {n_before} → {n_after}" +
          (f" ({dropped} pure-stop-codon exons dropped)" if dropped else ""))

    body = trimmed.select(GTF_COLUMNS).write_csv(separator="\t", include_header=False, quote_style="never")
    with open(out_path, "w") as fh:
        fh.write("##gff-version 2\n")
        fh.write("# Study 1 ORFs: VanHeesch et al. 2022 (PHASE I Ribo-seq ORFs)\n")
        fh.write("# Coordinates: 1-based (GTF format), GENCODE v45\n")
        fh.write("# CDS excludes stop codon (GENCODE convention)\n")
        fh.write(body)

    print(f"Study 1 GTF written to: {out_path}")


# ---------------------------------------------------------------------------
# Study 2
# ---------------------------------------------------------------------------

def process_study2(xlsx_path: Path, out_path: Path) -> None:
    df = pd.read_excel(xlsx_path, sheet_name="Supplementary Table 2", header=0)
    print(f"Study 2: loaded {len(df)} ORFs")
    df = df[df["primary_set"].astype(str).str.strip().str.lower() == "yes"].reset_index(drop=True)
    print(f"Study 2: {len(df)} ORFs after filtering to primary_set == 'yes'")

    rows: list[dict] = []
    for _, row in df.iterrows():
        chrom = f"chr{row['chrm']}"
        strand = str(row["strand"])

        # Comma-separated 0-based coordinates; convert to 1-based GTF
        starts_0based, ends_0based = parse_coords(
            row["starts (0-based)"], row["ends (0-based)"], sep=","
        )
        # BED: [start, end) 0-based  →  GTF: [start+1, end] 1-based
        starts_1based = [s + 1 for s in starts_0based]
        ends_1based = ends_0based  # BED end (exclusive, 0-based) == GTF end (inclusive, 1-based)

        rows.extend(build_orf_records(
            chrom=chrom, strand=strand,
            starts_1based=starts_1based, ends_1based=ends_1based,
            orf_id=str(row["releasev45_id"]),
            gene_id=str(row.get("gene_id", "")),
            gene_name=str(row.get("gene_name", "")),
            transcript_id=str(row.get("transcript", "")),
            orf_type=str(row.get("orf_type", "")),
            gene_biotype=str(row.get("gene_biotype", "")),
            source="study2_riboseq",
        ))

    gtf_df = pl.DataFrame(rows, schema={c: pl.String for c in GTF_COLUMNS})
    gtf_df = gtf_df.with_columns(
        pl.col("start").cast(pl.Int64),
        pl.col("end").cast(pl.Int64),
        pl.col("attributes").str.extract(r'orf_id "([^"]+)"').alias("transcript_id"),
    )

    n_before = gtf_df.filter(pl.col("feature") == "CDS").height
    trimmed = trim_stop_codon(gtf_df).drop("transcript_id")
    n_after = trimmed.filter(pl.col("feature") == "CDS").height
    dropped = n_before - n_after
    print(f"Study 2: CDS rows {n_before} → {n_after}" +
          (f" ({dropped} pure-stop-codon exons dropped)" if dropped else ""))

    body = trimmed.select(GTF_COLUMNS).write_csv(separator="\t", include_header=False, quote_style="never")
    with open(out_path, "w") as fh:
        fh.write("##gff-version 2\n")
        fh.write("# Study 2 ORFs: supplementary table 2\n")
        fh.write("# Coordinates: converted from 0-based BED to 1-based GTF, GENCODE v45\n")
        fh.write("# CDS excludes stop codon (GENCODE convention)\n")
        fh.write(body)

    print(f"Study 2 GTF written to: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    process_study2(
        xlsx_path=DATA_DIR / "media-1.xlsx",
        out_path=OUT_DIR / "study2_orfs.gtf",
    )
