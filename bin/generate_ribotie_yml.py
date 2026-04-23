#!/bin/env python
"""Generate a RiboTIE config YAML file from input file paths."""

import argparse
import csv
import os
import re
import yaml


BAM_SAMPLE_RE = re.compile(r'_([A-Za-z])_Unmapped')


def extract_sample_letter(bam_path):
    """Extract the single-letter sample identifier from a BAM filename."""
    basename = os.path.basename(bam_path)
    m = BAM_SAMPLE_RE.search(basename)
    if not m:
        raise ValueError(f"Could not extract sample letter from BAM filename: {basename}")
    return m.group(1)


def load_labels(csv_path):
    """Load sample labels CSV, returning a dict mapping letter -> (sample_name, condition)."""
    labels = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["Sample"].strip()
            condition = row["Condition"].strip()
            letter = sample.rsplit("_", 1)[-1]
            labels[letter] = (sample, condition)
    return labels


def main():
    parser = argparse.ArgumentParser(
        description="Generate a RiboTIE config YAML file."
    )
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--fa", required=True, help="Path to genome FASTA file")
    parser.add_argument(
        "--bams",
        nargs="+",
        required=True,
        help="Ribosome profiling BAM files",
    )
    parser.add_argument("--h5", required=True, help="Path to output HDF5 file")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--labels",
        required=True,
        help="CSV file with Sample and Condition columns",
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["merged", "separate"],
        help="'merged': all samples in one group; "
        "'separate': samples grouped by condition from --labels",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="ribotie_config.yml",
        help="Output YAML file path (default: ribotie_config.yml)",
    )

    args = parser.parse_args()

    label_map = load_labels(args.labels)

    ribo_paths = {}
    sample_condition = {}
    samples_by_condition = {}
    for bam in sorted(args.bams):
        letter = extract_sample_letter(bam)
        if letter not in label_map:
            raise ValueError(
                f"Letter '{letter}' from BAM '{os.path.basename(bam)}' not found in labels CSV"
            )
        sample_name, condition = label_map[letter]
        ribo_paths[sample_name] = bam
        sample_condition[sample_name] = condition
        samples_by_condition.setdefault(condition, []).append(sample_name)

    config = {
        "gtf_path": args.gtf,
        "fa_path": args.fa,
        "ribo_paths": ribo_paths,
        "h5_path": args.h5,
        "out_prefix": args.out_prefix,
        "samples": (
            {"merged": list(ribo_paths.keys())}
            if args.mode == "merged"
            else samples_by_condition
        ),
    }

    with open(args.output, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"Config written to {args.output}")
    print(f"Found {len(ribo_paths)} BAM files:")
    for sample, path in ribo_paths.items():
        print(f"  {sample} ({sample_condition[sample]}): {path}")


if __name__ == "__main__":
    main()
