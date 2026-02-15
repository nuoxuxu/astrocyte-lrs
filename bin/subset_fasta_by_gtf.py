#!/bin/env python
import argparse
import re


def get_transcript_ids_from_gtf(gtf_path):
    ids = set()
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m:
                ids.add(m.group(1))
    return ids


def subset_fasta(fasta_path, transcript_ids, output_path):
    with open(fasta_path) as fin, open(output_path, 'w') as fout:
        write = False
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                write = seq_id in transcript_ids
            if write:
                fout.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    transcript_ids = get_transcript_ids_from_gtf(args.gtf)
    subset_fasta(args.fasta, transcript_ids, args.output)
