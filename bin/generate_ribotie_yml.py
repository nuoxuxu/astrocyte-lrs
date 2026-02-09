#!/bin/env python
"""Generate a RiboTIE config YAML file from input file paths."""

import argparse
import glob
import os
import re
import yaml


BAM_SUFFIX = "_Unmapped.Aligned.toTranscriptome.out.bam"


def extract_sample_key(bam_path):
    """Extract sample key from BAM filename by removing the known suffix."""
    basename = os.path.basename(bam_path)
    if basename.endswith(BAM_SUFFIX):
        return basename[: -len(BAM_SUFFIX)]
    return os.path.splitext(basename)[0]


def parse_samples(sample_args, ribo_keys):
    """Parse sample groupings from 'GroupName:key1,key2' format."""
    samples = {}
    if not sample_args:
        return None
    for entry in sample_args:
        group_name, keys_str = entry.split(":", 1)
        keys = [k.strip() for k in keys_str.split(",")]
        for k in keys:
            if k not in ribo_keys:
                raise ValueError(
                    f"Sample key '{k}' not found in ribo_paths keys: {ribo_keys}"
                )
        samples[group_name] = keys
    return samples


def main():
    parser = argparse.ArgumentParser(
        description="Generate a RiboTIE config YAML file."
    )
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--fa", required=True, help="Path to genome FASTA file")
    parser.add_argument(
        "--bam-glob",
        required=True,
        help="Glob pattern for ribosome profiling BAM files "
        f"(files should end with '{BAM_SUFFIX}')",
    )
    parser.add_argument("--h5", required=True, help="Path to output HDF5 file")
    parser.add_argument("--out-prefix", required=True, help="Output prefix")
    parser.add_argument(
        "--samples",
        nargs="+",
        metavar="GROUP:key1,key2",
        help="Sample groupings, e.g. Unstim:astro_A,astro_B Stim:astro_C,astro_D",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="ribotie_config.yml",
        help="Output YAML file path (default: ribotie_config.yml)",
    )

    args = parser.parse_args()

    # Resolve BAM files from glob pattern
    bam_files = sorted(glob.glob(args.bam_glob))
    if not bam_files:
        parser.error(f"No BAM files matched glob pattern: {args.bam_glob}")

    # Build ribo_paths dict
    ribo_paths = {}
    for bam in bam_files:
        key = extract_sample_key(bam)
        ribo_paths[key] = bam

    # Build config
    config = {
        "gtf_path": args.gtf,
        "fa_path": args.fa,
        "ribo_paths": ribo_paths,
        "h5_path": args.h5,
        "out_prefix": args.out_prefix,
    }

    # Add samples if provided
    if args.samples:
        samples = parse_samples(args.samples, list(ribo_paths.keys()))
        config["samples"] = samples

    with open(args.output, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"Config written to {args.output}")
    print(f"Found {len(ribo_paths)} BAM files:")
    for key, path in ribo_paths.items():
        print(f"  {key}: {path}")


if __name__ == "__main__":
    main()