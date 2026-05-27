#!/usr/bin/env python
"""
Intersect novel coding splice junctions with BigBrain cis-sQTLs.

Steps:
  1. Load novel coding junctions from results/aim_2/novel_coding_junctions.tsv
  2. Parse from_collaborator/filtered_output.gtf to map each junction to the
     transcript_ids that contain it (semicolon-separated when multiple)
  3. Build feature key: {chrom}_{donor}_{acceptor-1}  (acceptor-1 converts the
     first base of the downstream exon to the last base of the intron, matching
     the BigBrain coordinate convention)
  4. Join on the feature column of BigBrain_cis_sQTL_ALL_top_assoc.tsv
  5. Keep feature, transcript_ids, gene/variant identifiers, and meta-analysis
     statistics → results/aim_2/novel_coding_junction_sqtl_matches.tsv
"""

import polars as pl

COLLAB_GTF     = "from_collaborator/filtered_output.gtf"
JUNCTIONS_FILE = "results/aim_2/novel_coding_junctions.tsv"
SQTL_FILE      = "data/BigBrain_cis_sQTL_ALL_top_assoc.tsv"
OUT_FILE       = "results/aim_2/novel_coding_junction_sqtl_matches.tsv"

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

    # Aggregate transcript_ids per junction
    return (
        junc_tx
        .group_by(["chrom", "donor", "acceptor", "strand"])
        .agg(pl.col("transcript_id").unique().sort().str.join(";").alias("transcript_ids"))
    )


print("Building junction → transcript_id map from collaborator GTF …")
junc_tx_map = build_junction_to_transcripts(COLLAB_GTF)

junctions = pl.read_csv(JUNCTIONS_FILE, separator="\t")

# Attach transcript_ids
junctions = junctions.join(junc_tx_map, on=["chrom", "donor", "acceptor", "strand"], how="left")

# Build feature key matching BigBrain convention: chrom_donor_(acceptor-1)
junctions = junctions.with_columns(
    (pl.col("chrom") + "_" + pl.col("donor").cast(pl.String) + "_" +
     (pl.col("acceptor") - 1).cast(pl.String)).alias("feature")
)

sqtl = pl.read_csv(SQTL_FILE, separator="\t").select(META_COLS)

matches = junctions.join(sqtl, on="feature", how="inner")

# Place transcript_ids right after feature for readability
col_order = ["chrom", "donor", "acceptor", "strand", "feature", "transcript_ids"] + \
            [c for c in matches.columns if c not in ("chrom", "donor", "acceptor", "strand", "feature", "transcript_ids")]
matches = matches.select(col_order)

print(f"Novel coding junctions:          {len(junctions):,}")
print(f"Matched to BigBrain sQTLs:       {len(matches):,}")

matches.write_csv(OUT_FILE, separator="\t")
print(f"Written to {OUT_FILE}")
