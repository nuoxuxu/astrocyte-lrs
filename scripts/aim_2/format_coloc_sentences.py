#!/usr/bin/env python3
"""
Generate one-sentence summaries for novel coding junction sQTL–GWAS colocalization hits.

For each unique (junction, disease) pair (keeping the highest PP4 row), prints:

  "sQTL [rsID] in gene [X] alters splicing of an astrocyte-reactive ncORF and
   colocalizes with [disease] GWAS (PP4 = [value]); the alternative splice event
   [creates an in-frame alteration in / introduces a frameshift in] the ncORF
   and is supported by [IF evidence]."

Frame effect is inferred from intron length (acceptor − donor − 1):
  divisible by 3 → in-frame alteration
  otherwise      → frameshift (NMD-prone)
"""

import csv
import sys
from pathlib import Path

COLOC_FILE = "nextflow_results/aim_2/novel_coding_junction_sqtl_coloc_matches.tsv"
MIN_PP4    = 0.5        # only report entries at or above this threshold


def first_gene(gene_names: str) -> str:
    """Return the first recognizable HGNC symbol from a semicolon-joined list."""
    for g in gene_names.split(";"):
        g = g.strip()
        if g and not g.startswith("ENSG"):
            return g
    return gene_names.split(";")[0].strip()


def frame_label(donor: int, acceptor: int) -> str:
    intron_length = acceptor - donor - 1
    if intron_length % 3 == 0:
        return "creates an in-frame alteration in"
    else:
        return "introduces a frameshift in"


def if_evidence(if_val: float) -> str:
    if if_val >= 0.5:
        return f"high isoform usage (IF = {if_val:.2f})"
    elif if_val >= 0.1:
        return f"moderate isoform usage (IF = {if_val:.2f})"
    elif if_val >= 0.01:
        return f"low isoform usage (IF = {if_val:.3f})"
    else:
        return f"very low isoform usage (IF = {if_val:.4f})"


def main():
    path = Path(COLOC_FILE)
    if not path.exists():
        sys.exit(f"File not found: {path}")

    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    # Deduplicate: keep highest PP4 per (feature, disease)
    best: dict[tuple, dict] = {}
    for row in rows:
        key = (row["feature"], row["disease"])
        pp4 = float(row["PP.H4.abf"])
        if key not in best or pp4 > float(best[key]["PP.H4.abf"]):
            best[key] = row

    hits = [r for r in best.values() if float(r["PP.H4.abf"]) >= MIN_PP4]
    hits.sort(key=lambda r: float(r["PP.H4.abf"]), reverse=True)

    print(f"Reporting {len(hits)} hits (PP4 >= {MIN_PP4})\n")

    for r in hits:
        gene    = first_gene(r["gene_names"])
        variant = r["variant_id"]
        disease = r["disease"]
        pp4     = float(r["PP.H4.abf"])
        donor   = int(r["donor"])
        acceptor = int(r["acceptor"])
        if_val  = float(r["IF"])

        splice = frame_label(donor, acceptor)
        evid   = if_evidence(if_val)

        sentence = (
            f"sQTL {variant} in gene {gene} alters splicing of an astrocyte-reactive "
            f"ncORF and colocalizes with {disease} GWAS (PP4 = {pp4:.3f}); "
            f"the alternative splice event {splice} the ncORF "
            f"and is supported by {evid}."
        )
        print(sentence)
        print()


if __name__ == "__main__":
    main()
