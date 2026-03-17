#!/usr/bin/env python3
"""
Generate a combined nf-core style samplesheet for ribo-seq and RNA-seq files.

Usage:
    python scripts/make_riboseq_samplesheet.py \
        --ribo_input data/ribo_seq \
        --rnaseq_input data/short_read \
        --rnaseq_primer data/short_read_primer.csv \
        --output samplesheet.csv
"""

import argparse
import csv
import re
from pathlib import Path

FIELDS = ["sample", "fastq_1", "fastq_2", "strandedness", "type"]


def riboseq_rows(input_dir: Path) -> list[dict]:
    rows = []
    for f in sorted(input_dir.iterdir()):
        if f.is_file():
            sample = re.sub(r"_Unmapped\.out\.mate1$", "", f.name)
            rows.append({
                "sample": sample,
                "fastq_1": str(f.resolve()),
                "fastq_2": "",
                "strandedness": "auto",
                "type": "riboseq",
            })
    return rows


def rnaseq_rows(input_dir: Path, primer_csv: Path) -> list[dict]:
    # Load prefix → sample mapping (deduplicate)
    prefix_to_sample: dict[str, str] = {}
    with primer_csv.open() as fh:
        for row in csv.DictReader(fh):
            prefix_to_sample[row["prefix"].strip()] = row["sample"].strip()

    # Group R1/R2 files by prefix
    pairs: dict[str, dict[str, Path]] = {}
    for f in sorted(input_dir.iterdir()):
        if not f.is_file():
            continue
        m = re.match(r"^(\S+?)_S\d+_L\d+_(R[12])_001\.fastq\.gz$", f.name)
        if not m:
            continue
        prefix, read = m.group(1), m.group(2)
        pairs.setdefault(prefix, {})[read] = f

    rows = []
    for prefix, reads in sorted(pairs.items()):
        sample = prefix_to_sample.get(prefix, prefix)
        rows.append({
            "sample": sample,
            "fastq_1": str(reads["R1"].resolve()) if "R1" in reads else "",
            "fastq_2": str(reads["R2"].resolve()) if "R2" in reads else "",
            "strandedness": "auto",
            "type": "rnaseq",
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description="Create combined nf-core samplesheet")
    parser.add_argument("--ribo_input", default="data/ribo_seq")
    parser.add_argument("--rnaseq_input", default="data/short_read")
    parser.add_argument("--rnaseq_primer", default="data/short_read_primer.csv")
    parser.add_argument("--output", default="samplesheet.csv")
    args = parser.parse_args()

    ribo_dir = Path(args.ribo_input)
    rna_dir = Path(args.rnaseq_input)
    primer_csv = Path(args.rnaseq_primer)

    for p in (ribo_dir, rna_dir, primer_csv):
        if not p.exists():
            raise SystemExit(f"Not found: {p}")

    rows = riboseq_rows(ribo_dir) + rnaseq_rows(rna_dir, primer_csv)

    output_path = Path(args.output)
    with output_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_path}")
    for r in rows:
        print(f"  [{r['type']:8s}] {r['sample']}")


if __name__ == "__main__":
    main()
