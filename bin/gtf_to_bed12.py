#!/usr/bin/env python3
"""Convert RiboTIE GTF (transcript/exon/CDS features) to BED12 format.

Produces one BED12 row per transcript, with:
  - exon blocks as BED12 blocks
  - CDS span as thickStart/thickEnd (renders coding region thick in UCSC browser)
  - ORF_id as the BED name field
  - ribotie_score (0-1) scaled to 0-1000 as BED score
"""

import sys
import re
from collections import defaultdict


def parse_attrs(attr_str):
    attrs = {}
    for m in re.finditer(r'(\w+)\s+"([^"]+)"', attr_str):
        attrs[m.group(1)] = m.group(2)
    return attrs


def main(gtf_file, bed_file):
    transcripts = {}

    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attr_str = parts
            start = int(start) - 1   # GTF is 1-based; BED is 0-based half-open
            end = int(end)
            attrs = parse_attrs(attr_str)
            tid = attrs.get('transcript_id', '')
            if not tid:
                continue

            if feature == 'transcript':
                t = transcripts.setdefault(tid, {
                    'chrom': chrom, 'start': start, 'end': end,
                    'strand': strand, 'exons': [], 'cds': [],
                    'orf_id': tid, 'orf_type': '.', 'ribotie_score': '0',
                    'gene_name': attrs.get('gene_name', '.'),
                })
                t['start'] = start
                t['end'] = end
                t['chrom'] = chrom
                t['strand'] = strand
                t['gene_name'] = attrs.get('gene_name', t.get('gene_name', '.'))

            elif feature == 'exon':
                t = transcripts.setdefault(tid, {
                    'chrom': chrom, 'start': start, 'end': end,
                    'strand': strand, 'exons': [], 'cds': [],
                    'orf_id': tid, 'orf_type': '.', 'ribotie_score': '0',
                    'gene_name': attrs.get('gene_name', '.'),
                })
                t['exons'].append((start, end))

            elif feature == 'CDS':
                t = transcripts.setdefault(tid, {
                    'chrom': chrom, 'start': start, 'end': end,
                    'strand': strand, 'exons': [], 'cds': [],
                    'orf_id': tid, 'orf_type': '.', 'ribotie_score': '0',
                    'gene_name': attrs.get('gene_name', '.'),
                })
                t['cds'].append((start, end))
                t['orf_id'] = attrs.get('ORF_id', tid)
                t['orf_type'] = attrs.get('ORF_type', '.')
                t['ribotie_score'] = attrs.get('ribotie_score', '0')

    with open(bed_file, 'w') as out:
        for tid, t in transcripts.items():
            exons = sorted(t['exons'])
            if not exons:
                continue

            tx_start = t['start']
            tx_end = t['end']

            if t['cds']:
                thick_start = min(c[0] for c in t['cds'])
                thick_end = max(c[1] for c in t['cds'])
            else:
                thick_start = tx_start
                thick_end = tx_start   # zero-width thick = UTR-only display

            try:
                bed_score = int(min(1000, max(0, float(t['ribotie_score']) * 1000)))
            except (ValueError, TypeError):
                bed_score = 0

            block_count = len(exons)
            block_sizes = ','.join(str(e - s) for s, e in exons) + ','
            block_starts = ','.join(str(s - tx_start) for s, e in exons) + ','

            out.write(
                f"{t['chrom']}\t{tx_start}\t{tx_end}\t{t['orf_id']}\t{bed_score}\t"
                f"{t['strand']}\t{thick_start}\t{thick_end}\t0\t"
                f"{block_count}\t{block_sizes}\t{block_starts}\n"
            )


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input.gtf> <output.bed>")
    main(sys.argv[1], sys.argv[2])
