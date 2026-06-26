#!/scratch/nxu/astrocyte-lrs/env/bin/python3
"""Add is_novel column to quality metrics.

is_novel = 1 if any CDS segment of the ORF overlaps a genomic region not
annotated as an exon in the GENCODE reference; 0 otherwise.

Usage:
    annotate_is_novel.py quality_metrics.tsv ribotie.gtf annotation.gtf -o output.tsv
"""
import argparse
import bisect
import csv
import re

import polars as pl


def parse_attrs(s):
    return {m.group(1): m.group(2) for m in re.finditer(r'(\w+) "([^"]*)"', s)}


def build_merged_exons(gtf_path):
    """Parse GTF exon features; return {chrom: sorted list of merged (start, end)}."""
    raw = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            raw.setdefault(fields[0], []).append((int(fields[3]), int(fields[4])))
    merged = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        result = [[ivs[0][0], ivs[0][1]]]
        for s, e in ivs[1:]:
            if s <= result[-1][1] + 1:
                result[-1][1] = max(result[-1][1], e)
            else:
                result.append([s, e])
        merged[chrom] = [(a, b) for a, b in result]
    return merged


def load_orf_cds(gtf_path):
    """Return {orf_id: (chrom, [(start, end), ...])} from CDS features in ribotie GTF."""
    segs = {}
    chroms = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = parse_attrs(fields[8])
            orf_id = attrs.get('ORF_id') or attrs.get('transcript_id')
            if not orf_id:
                continue
            chroms[orf_id] = fields[0]
            segs.setdefault(orf_id, []).append((int(fields[3]), int(fields[4])))
    return {oid: (chroms[oid], segs[oid]) for oid in segs}


def get_novel_regions(chrom, cds_segs, merged_exons, starts_cache):
    """Return list of (start, end) 1-based genomic coords in cds_segs not covered by merged exons."""
    intervals = merged_exons.get(chrom, [])
    starts = starts_cache.get(chrom, [])
    novel = []

    for seg_s, seg_e in cds_segs:
        idx = max(0, bisect.bisect_right(starts, seg_s) - 1)
        covered_to = seg_s - 1
        i = idx
        while i < len(intervals):
            iv_s, iv_e = intervals[i]
            if iv_s > seg_e:
                break
            if iv_e < seg_s:
                i += 1
                continue
            if iv_s > covered_to + 1:
                novel.append((max(covered_to + 1, seg_s), iv_s - 1))
            covered_to = max(covered_to, iv_e)
            i += 1

        if covered_to < seg_e:
            novel.append((max(covered_to + 1, seg_s), seg_e))

    return novel


def main():
    parser = argparse.ArgumentParser(
        description='Add is_novel column: 1 if ORF overlaps non-GENCODE-exon region.'
    )
    parser.add_argument('quality_metrics')
    parser.add_argument('ribotie_gtf')
    parser.add_argument('final_classification')
    parser.add_argument('annotation_gtf')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    merged_exons = build_merged_exons(args.annotation_gtf)
    starts_cache = {chrom: [iv[0] for iv in ivs] for chrom, ivs in merged_exons.items()}
    orf_cds = load_orf_cds(args.ribotie_gtf)

    clf = pl.read_parquet(args.final_classification).select(['isoform', 'structural_category'])
    struct_cat = dict(zip(clf['isoform'], clf['structural_category']))

    with open(args.quality_metrics) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t',
                                fieldnames=reader.fieldnames + ['structural_category', 'is_novel', 'novel_region', 'non_overlapping_dist'])
        writer.writeheader()
        for row in reader:
            orf_id = row['ORF_id']
            row['structural_category'] = struct_cat.get(row['transcript_id'], 'NA')
            if orf_id in orf_cds:
                chrom, segs = orf_cds[orf_id]
                regions = get_novel_regions(chrom, segs, merged_exons, starts_cache)
                row['is_novel'] = 1 if regions else 0
                row['novel_region'] = ';'.join(f'{chrom}:{s}-{e}' for s, e in regions) if regions else 'NA'
                row['non_overlapping_dist'] = sum(e - s + 1 for s, e in regions)
            else:
                row['is_novel'] = 'NA'
                row['novel_region'] = 'NA'
                row['non_overlapping_dist'] = 'NA'
            writer.writerow(row)


if __name__ == '__main__':
    main()
