#!/scratch/nxu/astrocyte-lrs/env/bin/python3
"""Add ribo_seq_coverage column: total Ribo-seq signal over novel regions / non_overlapping_dist.

For ORFs where is_novel=1, queries a Ribo-seq bigWig over the novel genomic regions
(already stored in the novel_region column) and divides the summed signal by
non_overlapping_dist. All other rows receive NA.

Usage:
    add_riboseq_coverage.py quality_metrics.tsv riboseq.bw \\
        --bigwig-tool /path/to/bigWigAverageOverBed -o output.tsv
"""
import argparse
import csv
import os
import subprocess
import tempfile


def parse_novel_regions(s):
    """Parse 'chr7:100-200;chr7:300-400' -> list of (chrom, start, end) 1-based inclusive."""
    regions = []
    for part in s.split(';'):
        chrom, coords = part.rsplit(':', 1)
        s_coord, e_coord = coords.split('-')
        regions.append((chrom, int(s_coord), int(e_coord)))
    return regions


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('quality_metrics')
    parser.add_argument('riboseq_bw')
    parser.add_argument('--bigwig-tool', required=True,
                        help='Path to bigWigAverageOverBed executable')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    # Read all rows; collect novel regions for BED construction
    rows = []
    novel_regions = {}   # orf_id -> [(chrom, s, e), ...]
    with open(args.quality_metrics) as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)
            if row['is_novel'] == '1':
                novel_regions[row['ORF_id']] = parse_novel_regions(row['novel_region'])

    # Write BED and run bigWigAverageOverBed
    bed_fd, bed_path = tempfile.mkstemp(suffix='.bed')
    bw_fd, bw_out_path = tempfile.mkstemp(suffix='.txt')
    os.close(bw_fd)
    try:
        with os.fdopen(bed_fd, 'w') as bed_f:
            for orf_id, regions in novel_regions.items():
                for i, (chrom, s, e) in enumerate(regions):
                    # bigWigAverageOverBed uses 0-based half-open coordinates
                    bed_f.write(f'{chrom}\t{s - 1}\t{e}\t{orf_id}__{i}\n')

        subprocess.run(
            [args.bigwig_tool, args.riboseq_bw, bed_path, bw_out_path],
            check=True, capture_output=True,
        )

        # Parse output: name, size, covered, sum, mean0, mean
        orf_signal = {}  # orf_id -> cumulative signal sum
        with open(bw_out_path) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                orf_id = parts[0].rsplit('__', 1)[0]
                try:
                    signal = float(parts[3])
                except (ValueError, IndexError):
                    signal = 0.0
                orf_signal[orf_id] = orf_signal.get(orf_id, 0.0) + signal
    finally:
        os.unlink(bed_path)
        os.unlink(bw_out_path)

    # Write output
    with open(args.output, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, delimiter='\t',
                                fieldnames=fieldnames + ['ribo_seq_coverage'])
        writer.writeheader()
        for row in rows:
            if row['is_novel'] == '1':
                dist = int(row['non_overlapping_dist'])
                total = orf_signal.get(row['ORF_id'], 0.0)
                row['ribo_seq_coverage'] = total / dist if dist > 0 else 'NA'
            else:
                row['ribo_seq_coverage'] = 'NA'
            writer.writerow(row)


if __name__ == '__main__':
    main()
