#!/usr/bin/env python
"""Label noncanonical ORFs as 'annotated' or 'novel' by comparing their protein
sequences against proteins translated from CDS regions of a reference GTF DataFrame.

A predicted ORF is 'annotated' if its protein sequence (stop codon stripped)
matches any protein translated from the reference CDS + genome FASTA. Otherwise
it is 'novel'.
"""
import sys
from collections import defaultdict

import polars as pl
import pyfaidx

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _translate(dna: str) -> str:
    """Translate DNA to protein, stopping at (and including) the first stop codon."""
    dna = dna.upper()
    aa = []
    for i in range(0, len(dna) - 2, 3):
        residue = CODON_TABLE.get(dna[i : i + 3], "X")
        aa.append(residue)
        if residue == "*":
            break
    return "".join(aa)


def _resolve_chrom(genome: pyfaidx.Fasta, seqname: str) -> str | None:
    """Return the FASTA key that matches *seqname*, trying with/without 'chr' prefix."""
    if seqname in genome:
        return seqname
    alt = seqname[3:] if seqname.startswith("chr") else f"chr{seqname}"
    if alt in genome:
        return alt
    return None


def extract_gtf_cds_proteins(gtf: pl.DataFrame, genome_fasta: str) -> dict[str, str]:
    """Translate CDS regions from a GTF DataFrame + *genome_fasta* into protein sequences.

    Args:
        gtf:          Polars DataFrame of a CDS-annotated GTF (e.g. loaded via gtfparse
                      or polars). Must have columns: seqname, feature, start, end, strand,
                      frame, and either a ``transcript_id`` column or an ``attributes``
                      column containing a ``transcript_id "..."`` field.
        genome_fasta: Path to the matching genome FASTA (must be indexed or indexable
                      by pyfaidx).

    Returns:
        Dict mapping transcript_id -> protein sequence (stop codon as '*' when present,
        omitted when the CDS ends without an in-frame stop, as is typical for GENCODE).
    """
    # Extract transcript_id from attributes string if not already a column
    if "transcript_id" not in gtf.columns:
        gtf = gtf.with_columns(
            pl.col("attributes").str.extract(r'transcript_id "([^"]+)"').alias("transcript_id")
        )

    gtf = (
        gtf
        .filter(pl.col("feature") == "CDS")
        .with_columns(
            pl.col("frame").cast(pl.Int32, strict=False).fill_null(0),
        )
        .filter(pl.col("transcript_id").is_not_null())
    )

    genome = pyfaidx.Fasta(genome_fasta)

    # Group CDS exons by transcript
    cds_by_tx: dict[str, list[dict]] = defaultdict(list)
    for row in gtf.iter_rows(named=True):
        cds_by_tx[row["transcript_id"]].append(row)

    proteins: dict[str, str] = {}
    skipped = 0
    for tx_id, exons in cds_by_tx.items():
        strand = exons[0]["strand"]
        chrom_key = _resolve_chrom(genome, exons[0]["seqname"])
        if chrom_key is None:
            skipped += 1
            continue

        # Sort in transcript 5'→3' order
        exons_sorted = sorted(exons, key=lambda x: x["start"], reverse=(strand == "-"))

        # Concatenate CDS sequence (GTF coords are 1-based inclusive)
        cds_seq = ""
        for exon in exons_sorted:
            chunk = str(genome[chrom_key][exon["start"] - 1 : exon["end"]])
            if strand == "-":
                chunk = _reverse_complement(chunk)
            cds_seq += chunk

        # Apply reading-frame offset reported for the first CDS exon
        frame = exons_sorted[0]["frame"] or 0
        proteins[tx_id] = _translate(cds_seq[frame:])

    if skipped:
        print(f"Skipped {skipped} transcripts: chromosome not found in FASTA", file=sys.stderr)

    return proteins


def label_noncanonical_orfs_by_sequence(
    gtf: pl.DataFrame,
    noncanonical_df: pl.DataFrame,
    genome_fasta: str,
    sequence_col: str = "orf_sequence",
) -> pl.DataFrame:
    """Label each row in *noncanonical_df* as 'annotated' or 'novel'.

    A row is 'annotated' if its protein sequence (stop codon stripped) matches
    any protein translated from the CDS regions of *gtf*. Otherwise it is
    'novel'.

    Args:
        gtf:             Polars DataFrame of a CDS-annotated reference GTF (e.g. GENCODE),
                         as loaded by gtfparse or polars.
        noncanonical_df: Polars DataFrame with at least:
                           - ``chrm``         chromosome (with or without 'chr' prefix)
                           - ``sequence_col`` protein sequence; stop codon denoted '*'
        genome_fasta:    Path to the genome FASTA.
        sequence_col:    Name of the column in *noncanonical_df* that contains protein
                         sequences. Default: ``'orf_sequence'``.

    Returns:
        *noncanonical_df* with an additional ``label`` column ('annotated' | 'novel').
    """
    ref_proteins = extract_gtf_cds_proteins(gtf, genome_fasta)

    # Build a set of reference sequences; strip trailing stop codon so the
    # comparison is robust to whether either side includes it
    ref_seqs: set[str] = {seq.rstrip("*") for seq in ref_proteins.values()}
    print(
        f"Reference: {len(ref_seqs)} unique protein sequences "
        f"from {len(ref_proteins)} transcripts",
        file=sys.stderr,
    )

    labels = [
        "annotated" if str(seq).rstrip("*") in ref_seqs else "novel"
        for seq in noncanonical_df[sequence_col]
    ]

    return noncanonical_df.with_columns(pl.Series("label", labels))
