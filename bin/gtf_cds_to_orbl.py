#!/usr/bin/env python3
"""Convert a GTF file with CDS features to the orbl.py input format.

Output columns: interval_string  strand  [extra_fields...]

where interval_string is chrN:start-end+chrN:start-end+... (sorted by start,
1-based inclusive coordinates as in GTF), and extra_fields are GTF attribute
values requested via --extra-fields.
"""

import sys
import argparse


DEFAULT_EXTRA_FIELDS = ['transcript_id', 'gene_name', 'gene_id', 'ribotie_score']


def parse_attributes(attr_string):
    attrs = {}
    for part in attr_string.strip().split(';'):
        part = part.strip()
        if not part:
            continue
        pieces = part.split('"')
        if len(pieces) >= 2:
            key = pieces[0].strip()
            val = pieces[1]
            attrs[key] = val
    return attrs


def main():
    parser = argparse.ArgumentParser(
        description='Convert GTF CDS features to orbl.py input format.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('gtf', help='Input GTF file')
    parser.add_argument(
        '--group-by',
        default='transcript_id',
        metavar='ATTR',
        help='GTF attribute to group CDS intervals by (default: transcript_id)',
    )
    parser.add_argument(
        '--extra-fields',
        nargs='+',
        default=DEFAULT_EXTRA_FIELDS,
        metavar='ATTR',
        help=(
            'GTF attributes to append as extra tab-separated fields after strand '
            '(default: %s)' % ' '.join(DEFAULT_EXTRA_FIELDS)
        ),
    )
    parser.add_argument(
        '--no-extra-fields',
        action='store_true',
        help='Output only interval string and strand, no extra fields',
    )
    args = parser.parse_args()

    extra_fields = [] if args.no_extra_fields else args.extra_fields

    # orf_data: group_key -> {'chrom': str, 'strand': str, 'intervals': set,
    #                          'attrs': dict}
    orf_data = {}
    insertion_order = []

    with open(args.gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attr_str = cols[:9]
            # Collect CDS intervals; also collect stop_codon intervals because GTF
            # CDS features exclude the stop codon but ORBL requires it as the last codon.
            if feature not in ('CDS', 'stop_codon'):
                continue

            attrs = parse_attributes(attr_str)
            group_key = attrs.get(args.group_by)
            if group_key is None:
                if feature == 'CDS':
                    print(
                        'WARNING: line missing attribute "%s", skipping: %s' % (args.group_by, line[:80]),
                        file=sys.stderr,
                    )
                continue

            istart, iend = int(start), int(end)

            if group_key not in orf_data:
                orf_data[group_key] = {
                    'chrom': chrom,
                    'strand': strand,
                    'intervals': set(),
                    'attrs': attrs,
                }
                insertion_order.append(group_key)

            orf_data[group_key]['intervals'].add((istart, iend))

    # Write output in GTF order
    for group_key in insertion_order:
        data = orf_data[group_key]
        # Sort intervals and merge adjacent/overlapping ones so the stop codon
        # is absorbed into the final CDS exon rather than listed separately.
        sorted_ivs = sorted(data['intervals'])
        merged = [list(sorted_ivs[0])]
        for s, e in sorted_ivs[1:]:
            if s <= merged[-1][1] + 1:  # adjacent or overlapping
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        interval_str = '+'.join(
            '%s:%d-%d' % (data['chrom'], s, e) for s, e in merged
        )
        out_cols = [interval_str, data['strand']]
        for field_name in extra_fields:
            out_cols.append(data['attrs'].get(field_name, '.'))
        print('\t'.join(out_cols))


if __name__ == '__main__':
    main()
