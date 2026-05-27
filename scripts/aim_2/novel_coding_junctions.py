#!/usr/bin/env python
"""
Identify novel coding splice junctions.

Steps:
  1. Parse CDS features from from_collaborator/filtered_output.gtf
  2. Parse CDS features from GENCODE v47 annotation GTF
  3. Derive intron junctions between consecutive CDS exons per transcript
  4. A junction is novel if the exact (chrom, donor, acceptor, strand) tuple
     is absent from GENCODE coding junctions. This includes junctions where
     both splice sites individually exist in GENCODE but are paired differently.
     → results/aim_2/novel_coding_junctions.tsv
"""

import polars as pl

COLLAB_GTF  = "from_collaborator/filtered_output.gtf"
GENCODE_GTF = "/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf"
OUT_FILE    = "results/aim_2/novel_coding_junctions.tsv"


def parse_attribute(attr_string: str, key: str) -> str:
    """Extract a single attribute value from a GTF attribute string."""
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

    # Deduplicate identical CDS entries (same transcript may have duplicate rows
    # with different ribotie_score or other per-read attributes)
    cds = cds.unique(subset=["chrom", "start", "end", "strand", "transcript_id"])

    # Sort by transcript then position (start, then end for stability)
    cds = cds.sort(["transcript_id", "start", "end"])

    # Derive junctions between consecutive CDS exons in the same transcript
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
        # Keep only true introns (gap between non-overlapping CDS exons)
        .filter(pl.col("donor") < pl.col("acceptor"))
        .unique()
    )
    return junctions


print("Parsing collaborator GTF …")
collab_junctions = load_cds_junctions(COLLAB_GTF)
print(f"  {len(collab_junctions):,} unique coding junctions")

print("Parsing GENCODE GTF …")
gencode_junctions = load_cds_junctions(GENCODE_GTF)
print(f"  {len(gencode_junctions):,} unique coding junctions")

# Set difference: in collaborator but not in GENCODE
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

novel.write_csv(OUT_FILE, separator="\t")
print(f"Written to {OUT_FILE}")
