#!/usr/bin/env python3
"""
Identify novel coding splice junctions.

Steps:
  1. Parse CDS features from a collaborator GTF
  2. Parse CDS features from GENCODE annotation GTF
  3. Derive intron junctions between consecutive CDS exons per transcript
  4. A junction is novel if the exact (chrom, donor, acceptor, strand) tuple
     is absent from GENCODE coding junctions.
"""

import argparse

import polars as pl


def parse_attribute(attr_string: str, key: str) -> str:
    for field in attr_string.split(";"):
        field = field.strip()
        if field.startswith(key + " "):
            return field[len(key) + 1:].strip().strip('"')
    return ""


def load_cds_junctions(gtf_path: str) -> pl.DataFrame:
    """
    Load CDS features from a GTF file and return a DataFrame of splice junctions
    (chrom, donor, acceptor, strand) derived from consecutive CDS exons per transcript.

    donor   = end coordinate of the upstream CDS exon (1-based inclusive)
    acceptor = start coordinate of the downstream CDS exon (1-based inclusive)
    """
    rows = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            chrom  = parts[0]
            start  = int(parts[3])
            end    = int(parts[4])
            strand = parts[6]
            attrs  = parts[8]
            tx_id  = parse_attribute(attrs, "transcript_id")
            rows.append((chrom, start, end, strand, tx_id))

    cds = pl.DataFrame(
        rows, schema=["chrom", "start", "end", "strand", "transcript_id"], orient="row"
    )

    cds = cds.unique(subset=["chrom", "start", "end", "strand", "transcript_id"])
    cds = cds.sort(["transcript_id", "start", "end"])

    shifted = cds.with_columns([
        pl.col("transcript_id").shift(-1).alias("next_tx"),
        pl.col("start").shift(-1).alias("next_start"),
    ])

    junctions = (
        shifted
        .filter(pl.col("transcript_id") == pl.col("next_tx"))
        .select([
            pl.col("chrom"),
            pl.col("end").alias("donor"),
            pl.col("next_start").alias("acceptor"),
            pl.col("strand"),
        ])
        .filter(pl.col("donor") < pl.col("acceptor"))
        .unique()
    )
    return junctions


def main():
    parser = argparse.ArgumentParser(description="Identify novel coding splice junctions.")
    parser.add_argument("predicted_cds_gtf", help="Collaborator GTF file")
    parser.add_argument("gencode_gtf", help="GENCODE annotation GTF file")
    parser.add_argument("-o", "--output", default="novel_coding_junctions.tsv",
                        help="Output TSV file (default: novel_coding_junctions.tsv)")
    args = parser.parse_args()

    print("Parsing collaborator GTF …")
    collab_junctions = load_cds_junctions(args.predicted_cds_gtf)
    print(f"  {len(collab_junctions):,} unique coding junctions")

    print("Parsing GENCODE GTF …")
    gencode_junctions = load_cds_junctions(args.gencode_gtf)
    print(f"  {len(gencode_junctions):,} unique coding junctions")

    gencode_set = set(
        zip(
            gencode_junctions["chrom"].to_list(),
            gencode_junctions["donor"].to_list(),
            gencode_junctions["acceptor"].to_list(),
            gencode_junctions["strand"].to_list(),
        )
    )

    novel = collab_junctions.filter(
        pl.struct(["chrom", "donor", "acceptor", "strand"]).map_elements(
            lambda r: (r["chrom"], r["donor"], r["acceptor"], r["strand"]) not in gencode_set,
            return_dtype=pl.Boolean,
        )
    )

    print(f"Novel coding junctions (collaborator − GENCODE): {len(novel):,}")

    novel.write_csv(args.output, separator="\t")
    print(f"Written to {args.output}")


if __name__ == "__main__":
    main()
