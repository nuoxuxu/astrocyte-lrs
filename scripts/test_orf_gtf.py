#!/usr/bin/env python3
"""
Validate the generated ORF GTF files by:
  1. Parsing the CDS blocks from the GTF
  2. Extracting genomic sequences with pyfaidx
  3. Translating to protein
  4. Comparing with reference protein sequences from the original Excel files

Run from the repo root:
  conda run -n patch_seq_spl python3 scripts/test_orf_gtf.py
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pyfaidx

GENOME_FA = Path("/Users/xunuo/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa")
DATA_DIR = Path("data")

# Standard codon table (DNA)
CODON_TABLE: dict[str, str] = {
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

RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def translate(nt_seq: str) -> str:
    """Translate a nucleotide sequence to protein (including stop codon as *)."""
    protein = []
    for i in range(0, len(nt_seq) - 2, 3):
        codon = nt_seq[i : i + 3].upper()
        aa = CODON_TABLE.get(codon, "?")
        protein.append(aa)
        if aa == "*":
            break
    return "".join(protein)


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf_cds(gtf_path: Path) -> dict[str, dict]:
    """
    Parse CDS records from a GTF file.

    Returns a dict keyed by orf_id:
      {
        "chrom": str,
        "strand": str,
        "cds_blocks": [(start_1based, end_1based), ...],  # sorted ascending
      }
    """
    orf_id_re = re.compile(r'orf_id "([^"]+)"')
    records: dict[str, dict] = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue

            chrom, _, _, start, end, _, strand, _, attrs = fields
            m = orf_id_re.search(attrs)
            if not m:
                continue
            orf_id = m.group(1)

            if orf_id not in records:
                records[orf_id] = {"chrom": chrom, "strand": strand, "cds_blocks": []}
            records[orf_id]["cds_blocks"].append((int(start), int(end)))

    # Sort blocks ascending by start (we handle strand during extraction)
    for rec in records.values():
        rec["cds_blocks"].sort(key=lambda x: x[0])

    return records


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def extract_orf_sequence(
    fa: pyfaidx.Fasta,
    chrom: str,
    strand: str,
    cds_blocks: list[tuple[int, int]],
) -> str:
    """
    Extract and concatenate CDS blocks from the genome.
    Blocks should be sorted ascending. pyfaidx uses 0-based half-open intervals.
    """
    # pyfaidx: fa[chrom][start:end] with 0-based half-open coords
    # GTF: 1-based inclusive  →  pyfaidx: start-1, end
    if strand == "+":
        parts = [str(fa[chrom][start - 1 : end]) for start, end in cds_blocks]
        return "".join(parts)
    else:
        # Process blocks right to left (3'→5' on genome = 5'→3' on mRNA)
        parts = [
            reverse_complement(str(fa[chrom][start - 1 : end]))
            for start, end in reversed(cds_blocks)
        ]
        return "".join(parts)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def normalize_aa(seq: str) -> str:
    """Strip trailing stop codon for comparison."""
    return seq.rstrip("*") if isinstance(seq, str) else ""


def validate_study(
    gtf_path: Path,
    gtf_records: dict[str, dict],
    ref_df: pd.DataFrame,
    orf_id_col: str,
    aa_col: str,
    fa: pyfaidx.Fasta,
    study_label: str,
) -> None:
    print(f"\n{'='*60}")
    print(f"  {study_label}")
    print(f"{'='*60}")
    print(f"  GTF: {gtf_path}")
    print(f"  ORFs in GTF: {len(gtf_records)}")

    # Build reference dict: orf_id -> aa_sequence
    ref_df = ref_df.dropna(subset=[aa_col])
    ref = dict(zip(ref_df[orf_id_col].astype(str), ref_df[aa_col].astype(str)))
    print(f"  Reference ORFs with AA sequence: {len(ref)}")

    STOP_CODONS = {"TAA", "TAG", "TGA"}
    total = match = mismatch = skipped = stop_at_end = 0

    mismatches: list[dict] = []

    for orf_id, rec in gtf_records.items():
        ref_aa = ref.get(orf_id)
        if ref_aa is None:
            skipped += 1
            continue

        total += 1
        try:
            nt_seq = extract_orf_sequence(fa, rec["chrom"], rec["strand"], rec["cds_blocks"])
            pred_aa = translate(nt_seq)
        except Exception as e:
            print(f"  ERROR extracting {orf_id}: {e}")
            skipped += 1
            continue

        last_codon = nt_seq[-3:].upper() if len(nt_seq) >= 3 else ""
        if last_codon in STOP_CODONS:
            stop_at_end += 1

        if normalize_aa(pred_aa) == normalize_aa(ref_aa):
            match += 1
        else:
            mismatch += 1
            if len(mismatches) < 10:  # keep first 10 examples
                mismatches.append(
                    {
                        "orf_id": orf_id,
                        "chrom": rec["chrom"],
                        "strand": rec["strand"],
                        "cds_blocks": rec["cds_blocks"],
                        "predicted": pred_aa,
                        "reference": ref_aa,
                    }
                )

    print(f"\n  Results ({total} ORFs with reference):")
    print(f"    Match:    {match:6d}  ({100*match/total:.1f}%)")
    print(f"    Mismatch: {mismatch:6d}  ({100*mismatch/total:.1f}%)")
    print(f"    Stop codon at end of CDS: {stop_at_end}/{total}"
          + (" ← EXPECTED 0 after trimming" if stop_at_end == 0 else " ← WARNING: stop codons still present"))
    print(f"    Skipped (no reference AA): {skipped}")

    if mismatches:
        print(f"\n  First {len(mismatches)} mismatch example(s):")
        for ex in mismatches:
            print(f"\n    orf_id  : {ex['orf_id']}")
            print(f"    location: {ex['chrom']}:{ex['strand']} blocks={ex['cds_blocks']}")
            print(f"    predicted : {ex['predicted'][:80]}")
            print(f"    reference : {ex['reference'][:80]}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Loading genome: {GENOME_FA}")
    fa = pyfaidx.Fasta(str(GENOME_FA))

    # --- Study 1 ---
    gtf1 = DATA_DIR / "study1_orfs.gtf"
    records1 = parse_gtf_cds(gtf1)
    df1 = pd.read_excel(
        DATA_DIR / "41587_2022_1369_MOESM2_ESM.xlsx",
        sheet_name="S2. PHASE I Ribo-seq ORFs",
        header=0,
    )
    validate_study(
        gtf_path=gtf1,
        gtf_records=records1,
        ref_df=df1,
        orf_id_col="orf_name",
        aa_col="orf_sequence",
        fa=fa,
        study_label="Study 1 – VanHeesch et al. 2022 (PHASE I Ribo-seq ORFs)",
    )

    # --- Study 2 ---
    gtf2 = DATA_DIR / "study2_orfs.gtf"
    records2 = parse_gtf_cds(gtf2)
    df2 = pd.read_excel(
        DATA_DIR / "media-1.xlsx",
        sheet_name="Supplementary Table 2",
        header=0,
    )
    # Use sequence_aa as primary reference; fall back to sequence_aa_MS if present
    aa_col = "sequence_aa"
    validate_study(
        gtf_path=gtf2,
        gtf_records=records2,
        ref_df=df2,
        orf_id_col="releasev45_id",
        aa_col=aa_col,
        fa=fa,
        study_label="Study 2 – Supplementary Table 2",
    )


if __name__ == "__main__":
    main()
