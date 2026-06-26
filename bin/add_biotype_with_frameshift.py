#!/usr/bin/env python3
"""Add biotype_with_frameshift column to quality_metrics.tsv for ORBLq calculation.

Derives biotype_with_frameshift from ORF_type_ORFanage and the pre-computed
frame_wrt_canonical_TIS column in the RiboTIE CSV:

    frame = frame_wrt_canonical_TIS  (0, 1, or 2)

Overlapping types (uoORF, intORF, doORF) are annotated as e.g. uoORF+1 or uoORF+2;
frame 0 (in-frame with canonical CDS) yields the bare type name without a suffix.
Types with no frameshift variant (uORF, dORF, lncRNA-ORF) are passed through directly.
All other types get NA.

BiotypesWithFS = ['uORF', 'uoORF+1', 'uoORF+2', 'intORF+1', 'intORF+2',
                  'doORF+1', 'doORF+2', 'dORF', 'lncRNA-ORF']

Usage:
    add_biotype_with_frameshift.py quality_metrics.tsv ribotie.csv -o output.tsv
"""
import argparse
import csv

DIRECT = {'uORF': 'uORF', 'dORF': 'dORF', 'lncRNA-ORF': 'lncRNA-ORF'}
NEEDS_FRAME = {'uoORF', 'intORF', 'doORF'}


def load_frame_map(ribotie_csv):
    """Return ORF_id -> frame (int 0/1/2 or None) from pre-computed ribotie column."""
    frame_map = {}
    with open(ribotie_csv) as f:
        for row in csv.DictReader(f):
            orf_id = row['ORF_id']
            raw = row.get('frame_wrt_canonical_TIS', '')
            try:
                frame_map[orf_id] = int(float(raw))
            except (ValueError, TypeError):
                frame_map[orf_id] = None
    return frame_map


def main():
    parser = argparse.ArgumentParser(
        description='Add biotype_with_frameshift to quality_metrics.tsv.'
    )
    parser.add_argument('quality_metrics')
    parser.add_argument('ribotie_csv')
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    frame_map = load_frame_map(args.ribotie_csv)

    with open(args.quality_metrics) as f_in, open(args.output, 'w', newline='') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        writer = csv.DictWriter(f_out, delimiter='\t',
                                fieldnames=reader.fieldnames + ['biotype_with_frameshift'])
        writer.writeheader()
        for row in reader:
            orf_id = row['ORF_id']
            orf_type = row['ORF_type_ORFanage']

            if orf_type in DIRECT:
                biotype_fs = DIRECT[orf_type]
            elif orf_type in NEEDS_FRAME:
                frame = frame_map.get(orf_id)
                biotype_fs = f'{orf_type}+{frame}' if frame in (1, 2) else orf_type
            else:
                biotype_fs = 'NA'

            row['biotype_with_frameshift'] = biotype_fs
            writer.writerow(row)


if __name__ == '__main__':
    main()
