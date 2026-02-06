#!/usr/bin/env python3
"""
Convert chromosome names in a chrom.sizes file to match those used in bedGraph files
using the mapping from hg38.chromAlias.txt
"""

import sys
import argparse
import glob

def parse_chrom_alias(alias_file):
    """
    Parse chromAlias file to create mappings between UCSC names and aliases.
    """
    alias_to_ucsc = {}
    ucsc_to_aliases = {}
    
    with open(alias_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            ucsc_name = parts[0]

            # Map all aliases (including the UCSC name itself) to the UCSC name
            for alias in parts:
                alias_to_ucsc[alias] = ucsc_name
                ucsc_to_aliases.setdefault(ucsc_name, set()).add(alias)

    return alias_to_ucsc, ucsc_to_aliases

def get_bedgraph_contigs(bedgraph_glob):
    """
    Extract unique contig names from all bedGraph files matching glob pattern
    """
    contigs = set()
    files = glob.glob(bedgraph_glob)
    
    if not files:
        print(f"Warning: No files match pattern '{bedgraph_glob}'", file=sys.stderr)
        return contigs
    
    print(f"Found {len(files)} files matching pattern", file=sys.stderr)
    
    for bedgraph_file in files:
        print(f"  Processing: {bedgraph_file}", file=sys.stderr)
        with open(bedgraph_file, 'r') as f:
            for line in f:
                if line.startswith('track') or line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    contigs.add(parts[0])
    
    return contigs

def convert_chrom_sizes(chrom_sizes_file, alias_to_ucsc, bedgraph_contigs, output_file):
    """
    Read chrom.sizes and output a new version with contig names converted
    to those used in the bedGraph files.
    Only include chromosomes that are present in the bedGraph files.
    """
    # Map bedGraph contigs to UCSC names
    bedgraph_by_ucsc = {}
    unmapped_contigs = []

    for contig in bedgraph_contigs:
        if contig in alias_to_ucsc:
            ucsc_name = alias_to_ucsc[contig]
        else:
            # If not found in alias file, keep it as is (might already be UCSC name)
            ucsc_name = contig
            unmapped_contigs.append(contig)

        bedgraph_by_ucsc.setdefault(ucsc_name, set()).add(contig)

    if unmapped_contigs:
        print(f"Warning: {len(unmapped_contigs)} contigs not found in alias file:", file=sys.stderr)
        for contig in sorted(unmapped_contigs)[:10]:
            print(f"  {contig}", file=sys.stderr)
        if len(unmapped_contigs) > 10:
            print(f"  ... and {len(unmapped_contigs) - 10} more", file=sys.stderr)

    # Read chrom.sizes and output only matching chromosomes, using bedGraph contig names
    written_count = 0

    with open(chrom_sizes_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_name = parts[0]

                if chrom_name in bedgraph_by_ucsc:
                    contig_names = sorted(bedgraph_by_ucsc[chrom_name])
                    # Prefer UCSC name if it exists in bedGraph; otherwise use first matching alias
                    if chrom_name in bedgraph_by_ucsc[chrom_name]:
                        out_name = chrom_name
                    else:
                        out_name = contig_names[0]

                    f_out.write(f"{out_name}\t{parts[1]}\n")
                    written_count += 1

    print(f"Converted chrom.sizes written to {output_file}", file=sys.stderr)
    print(f"Included {written_count} chromosomes out of {len(bedgraph_contigs)} bedGraph contigs", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description='Convert chrom.sizes to match bedGraph contig names using chromAlias mapping'
    )
    parser.add_argument('bedgraph_pattern', help='Glob pattern for bedGraph files (e.g., "*.bedgraph")')
    parser.add_argument('chrom_sizes', help='Input chrom.sizes file')
    parser.add_argument('chrom_alias', help='Chromosome alias mapping file (hg38.chromAlias.txt)')
    parser.add_argument('output', help='Output chrom.sizes file')
    
    args = parser.parse_args()
    
    # Parse the alias mapping
    print("Parsing chromosome alias file...", file=sys.stderr)
    alias_to_ucsc, _ = parse_chrom_alias(args.chrom_alias)
    print(f"Loaded {len(alias_to_ucsc)} alias mappings", file=sys.stderr)
    
    # Get contigs from all bedGraph files matching pattern
    print(f"Extracting contigs from bedGraph files matching '{args.bedgraph_pattern}'...", file=sys.stderr)
    bedgraph_contigs = get_bedgraph_contigs(args.bedgraph_pattern)
    print(f"Found {len(bedgraph_contigs)} unique contigs in bedGraph files", file=sys.stderr)
    
    # Convert chrom.sizes
    print("Converting chrom.sizes...", file=sys.stderr)
    convert_chrom_sizes(args.chrom_sizes, alias_to_ucsc, bedgraph_contigs, args.output)

if __name__ == '__main__':
    main()
