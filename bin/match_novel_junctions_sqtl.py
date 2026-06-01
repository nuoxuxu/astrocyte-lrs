#!/usr/bin/env python3
"""
Intersect novel coding splice junctions with BigBrain cis-sQTLs.

Steps:
  1. Load novel coding junctions TSV
  2. Parse predicted CDS GTF to map each junction to the transcript_ids that
     contain it (semicolon-separated when multiple)
  3. Build feature key: {chrom}_{donor}_{acceptor-1}  (acceptor-1 converts the
     first base of the downstream exon to the last base of the intron, matching
     the BigBrain coordinate convention)
  4. Join on the feature column of the BigBrain cis-sQTL file
  5. Add IF column: max IF_overall across all transcripts carrying the junction
  6. Keep feature, transcript_ids, IF, gene/variant identifiers, and meta-analysis
     statistics
"""

import argparse

import polars as pl


META_COLS = [
    "feature", "gene_ids", "gene_names",
    "variant_id", "chr", "pos", "ref", "alt", "Allele",
    "fixed_beta", "fixed_sd", "fixed_z", "Random_Z",
    "Fixed_P", "Random_P", "Fixed_bonf", "Random_bonf",
    "Fixed_FDR", "Random_FDR", "Bonf_FDR",
]


def parse_attribute(attr_string: str, key: str) -> str:
    for field in attr_string.split(";"):
        field = field.strip()
        if field.startswith(key + " "):
            return field[len(key) + 1:].strip().strip('"')
    return ""


def build_junction_to_transcripts(gtf_path: str) -> pl.DataFrame:
    """Return a DataFrame mapping (chrom, donor, acceptor, strand) → transcript_ids."""
    rows = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            rows.append((
                parts[0],
                int(parts[3]),
                int(parts[4]),
                parts[6],
                parse_attribute(parts[8], "transcript_id"),
            ))

    cds = pl.DataFrame(
        rows, schema=["chrom", "start", "end", "strand", "transcript_id"], orient="row"
    )
    cds = cds.unique(subset=["chrom", "start", "end", "strand", "transcript_id"])
    cds = cds.sort(["transcript_id", "start", "end"])

    shifted = cds.with_columns([
        pl.col("transcript_id").shift(-1).alias("next_tx"),
        pl.col("start").shift(-1).alias("next_start"),
    ])

    junc_tx = (
        shifted
        .filter(pl.col("transcript_id") == pl.col("next_tx"))
        .select([
            pl.col("chrom"),
            pl.col("end").alias("donor"),
            pl.col("next_start").alias("acceptor"),
            pl.col("strand"),
            pl.col("transcript_id"),
        ])
        .filter(pl.col("donor") < pl.col("acceptor"))
    )

    return (
        junc_tx
        .group_by(["chrom", "donor", "acceptor", "strand"])
        .agg(pl.col("transcript_id").unique().sort().str.join(";").alias("transcript_ids"))
    )


def main():
    parser = argparse.ArgumentParser(
        description="Intersect novel coding splice junctions with BigBrain cis-sQTLs."
    )
    parser.add_argument("predicted_cds_gtf", help="Predicted CDS GTF file")
    parser.add_argument("junctions", help="Novel coding junctions TSV (from novel_coding_junctions.py)")
    parser.add_argument("sqtl_file", help="BigBrain cis-sQTL TSV file")
    parser.add_argument("iso_file", help="IsoformSwitchAnalyzeR isoformFeatures.csv")
    parser.add_argument("-o", "--output", default="novel_coding_junction_sqtl_matches.tsv",
                        help="Output TSV file (default: novel_coding_junction_sqtl_matches.tsv)")
    args = parser.parse_args()

    print("Building junction → transcript_id map from predicted CDS GTF …")
    junc_tx_map = build_junction_to_transcripts(args.predicted_cds_gtf)

    junctions = pl.read_csv(args.junctions, separator="\t")

    junctions = junctions.join(junc_tx_map, on=["chrom", "donor", "acceptor", "strand"], how="left")

    junctions = junctions.with_columns(
        (pl.col("chrom") + "_" + pl.col("donor").cast(pl.String) + "_" +
         (pl.col("acceptor") - 1).cast(pl.String)).alias("feature")
    )

    sqtl = pl.read_csv(args.sqtl_file, separator="\t").select(META_COLS)

    matches = junctions.join(sqtl, on="feature", how="inner")

    iso = pl.read_csv(args.iso_file, null_values="NA")
    iso_map = dict(zip(
        iso["isoform_id"].to_list(),
        iso["IF_overall"].to_list(),
    ))

    def max_if(tx_ids: str) -> float | None:
        vals = [iso_map.get(tx.strip()) for tx in tx_ids.split(";")]
        vals = [v for v in vals if v is not None]
        return max(vals) if vals else None

    matches = matches.with_columns(
        pl.col("transcript_ids").map_elements(max_if, return_dtype=pl.Float64).alias("IF")
    )

    col_order = ["chrom", "donor", "acceptor", "strand", "feature", "transcript_ids", "IF"] + \
                [c for c in matches.columns if c not in ("chrom", "donor", "acceptor", "strand", "feature", "transcript_ids", "IF")]
    matches = matches.select(col_order)

    print(f"Novel coding junctions:          {len(junctions):,}")
    print(f"Matched to BigBrain sQTLs:       {len(matches):,}")

    matches.write_csv(args.output, separator="\t")
    print(f"Written to {args.output}")


if __name__ == "__main__":
    main()
