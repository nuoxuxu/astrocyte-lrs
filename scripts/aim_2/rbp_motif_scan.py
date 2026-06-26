#!/usr/bin/env python
"""
Scan sequences flanking sQTL variants for disruption of RBP binding motifs.

Motifs are IUPAC RNA consensus sequences drawn from:
  - EuPRI / CisBP-RNA (Sasse et al. 2024/2025, Nature Biotechnology)
  - Ray et al. 2013 (Nature) – original RNAcompete compendium
  - van Nostrand et al. 2020 (Nature) – ENCODE eCLIP consensus motifs

Approach:
  For each variant, extract ±40 bp genomic sequence (REF allele),
  build ALT sequence by swapping the central base, convert both to RNA,
  then find all IUPAC motif matches in each. Report motifs that are
  present in REF but absent (or weaker) in ALT (disrupted) or absent
  in REF but present in ALT (created).
"""

import re
import subprocess

GENOME = "/scratch/nxu/astrocytes/data/GRCh38.primary_assembly.genome.fa"
SAMTOOLS = "/scratch/nxu/astrocyte-lrs/env/bin/samtools"
FLANK = 40  # bp on each side of the variant

# ── IUPAC helpers ────────────────────────────────────────────────────────────
IUPAC = {
    "A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
    "R": "[AG]", "Y": "[CU]", "S": "[CG]", "W": "[AU]",
    "M": "[AC]", "K": "[GU]", "H": "[ACU]", "B": "[CGU]",
    "V": "[ACG]", "D": "[AGU]", "N": "[ACGU]",
}

def iupac_to_regex(motif):
    """Convert IUPAC RNA motif to regex pattern (case-insensitive on RNA)."""
    return "".join(IUPAC.get(b.upper(), b) for b in motif)

def dna_to_rna(seq):
    return seq.upper().replace("T", "U")

def find_motif_hits(rna_seq, motif_regex, motif_len):
    """Return start positions (0-based) of all non-overlapping matches."""
    return [m.start() for m in re.finditer(f"(?=({motif_regex}))", rna_seq)]

# ── RBP motif catalogue (curated from EuPRI / CisBP-RNA / ENCODE eCLIP) ────
# Format: (RBP_name, IUPAC_RNA_motif, role_in_splicing)
MOTIFS = [
    # SR proteins – exonic / intronic splicing enhancers
    ("SRSF1",    "GAAGAAGGA",  "ESE"),
    ("SRSF1",    "RGAAGAAC",   "ESE"),
    ("SRSF2",    "GCGYGG",     "ESE"),
    ("SRSF3",    "CWWCWY",     "ESE"),
    ("SRSF5",    "GGAUGAAG",   "ESE"),
    ("SRSF6",    "GAAGGA",     "ESE"),
    ("SRSF7",    "GACGAYGAY",  "ESE"),
    # hnRNPs – mostly silencers
    ("HNRNPA1",  "UAGGGU",     "ISS/ESS"),
    ("HNRNPA2",  "UAGG",       "ISS/ESS"),
    ("HNRNPC",   "UUUUU",      "ISS"),
    ("HNRNPF",   "AGGGGA",     "ISS"),
    ("HNRNPH1",  "GGGGA",      "ISS"),
    ("HNRNPH2",  "GGGGA",      "ISS"),
    ("HNRNPK",   "CCCC",       "ISS"),
    ("HNRNPL",   "CACAAG",     "ISE/ISS"),
    ("HNRNPM",   "UGUGUGUG",   "ISS"),
    # Neuronally-enriched / brain splicing factors
    ("RBFOX1",   "UGCAUG",     "ISE/ESE"),
    ("RBFOX2",   "UGCAUG",     "ISE/ESE"),
    ("NOVA1",    "UCAUY",      "ISE/ISS"),
    ("NOVA2",    "UCAUY",      "ISE/ISS"),
    ("PTBP1",    "UCUUC",      "ISS/ESS"),
    ("PTBP1",    "UCUCU",      "ISS/ESS"),
    ("PTBP2",    "UCUUC",      "ISS/ESS"),
    ("MBNL1",    "YGCY",       "ISE/ESE"),
    ("MBNL2",    "YGCY",       "ISE/ESE"),
    ("QKI",      "ACUAAY",     "ISE"),
    ("MSI1",     "RNAHM",      "ISE"),   # (N=any, H=non-G, M=AC)
    # ALS/FTD-relevant RBPs
    ("TARDBP",   "GUGUGUG",    "ISS"),
    ("TARDBP",   "UGUGUGU",    "ISS"),
    ("FUS",      "GUGGU",      "ISE"),
    ("EWSR1",    "UGGGGU",     "ISE"),
    # General post-transcriptional regulators
    ("ELAVL1",   "UUUUUU",     "stability/splicing"),
    ("ELAVL1",   "AUUUA",      "stability"),
    ("ELAVL2",   "UUUUUU",     "stability/splicing"),
    ("CELF1",    "UGUUUGU",    "ISS/ESS"),
    ("CELF2",    "UGUUUGU",    "ISS/ESS"),
    ("ZC3H12A",  "AUUUA",      "mRNA decay"),
    # Splice-site proximal factors
    ("U2AF2",    "UUUUUU",     "3'ss polypyrimidine"),
    ("U2AF1",    "YAG",        "3'ss AG"),
    ("SF1",      "YURAY",      "branch point"),
    # m6A writers/readers (affect splicing)
    ("YTHDF1",   "GGACU",      "m6A reader"),
    ("YTHDC1",   "GAACU",      "m6A-splicing"),
]

# ── Variant definitions ───────────────────────────────────────────────────────
VARIANTS = [
    # (label, gene, disease, chrom, pos_1based, ref_allele, alt_allele)
    # Primary: leafcutter sQTL top variant (pos from GRCh38-aligned analysis)
    # DGKZ×AD: rs57591833 (top sQTL, PP.H4=0.753); rs10838542 (second intron, PP.H4=0.518)
    ("DGKZ×AD (rs57591833)",    "DGKZ",    "AD",  "chr11", 45592862, "A", "G"),
    ("DGKZ×AD (rs10838542)",    "DGKZ",    "AD",  "chr11", 46023041, "A", "G"),
    # EFEMP1×SCZ: rs60772612 (PP.H4=0.708)
    ("EFEMP1×SCZ (rs60772612)", "EFEMP1",  "SCZ", "chr2",  55882632, "A", "T"),
    # SEPTIN4×AD: rs9908667 (PP.H4=0.537)
    ("SEPTIN4×AD (rs9908667)",  "SEPTIN4", "AD",  "chr17", 58300106, "C", "G"),
]

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

def rc(seq):
    return seq.translate(COMPLEMENT)[::-1]

def get_sequence(chrom, pos_1based, flank):
    """Fetch genomic sequence via samtools faidx (1-based coordinates)."""
    region = f"{chrom}:{pos_1based - flank}-{pos_1based + flank}"
    result = subprocess.run(
        [SAMTOOLS, "faidx", GENOME, region],
        capture_output=True, text=True, check=True
    )
    lines = result.stdout.strip().split("\n")[1:]  # skip header
    return "".join(lines).upper()

# ── Build regex patterns once ─────────────────────────────────────────────────
COMPILED = [(rbp, iupac_to_regex(motif), motif, role)
            for rbp, motif, role in MOTIFS]

def scan(rna_seq, snp_idx):
    """Return list of (rbp, motif, role, positions) that overlap snp_idx."""
    hits = []
    for rbp, regex, iupac, role in COMPILED:
        mlen = len(iupac.replace("(","").replace(")",""))
        for m in re.finditer(f"(?=({regex}))", rna_seq):
            start = m.start()
            end   = start + mlen
            if start <= snp_idx < end:
                hits.append((rbp, iupac, role, start, end))
    return hits

# ── Main ──────────────────────────────────────────────────────────────────────
results = {}

for label, gene, disease, chrom, pos, ref, alt in VARIANTS:
    dna_ref_plus = get_sequence(chrom, pos, FLANK)
    snp_idx = FLANK  # variant is at the centre of the window

    # BigBrain file uses "ref/alt" as regression alleles, not VCF convention.
    # GRCh38 may carry either allele. We try, in order:
    #   1. Plus strand, file-ref → file-alt substitution
    #   2. Plus strand, file-alt → file-ref substitution (GRCh38 carries alt)
    #   3. Minus strand (reverse complement), same logic
    genomic_base = dna_ref_plus[snp_idx]
    dna_ref_rc = rc(dna_ref_plus)
    genomic_base_minus = dna_ref_rc[snp_idx]

    if genomic_base == ref:
        dna_ref = dna_ref_plus
        dna_alt = dna_ref[:snp_idx] + alt + dna_ref[snp_idx+1:]
        strand = "+"
    elif genomic_base == alt:
        # GRCh38 carries the alt allele; swap so ref=genome, subst to reported-ref
        dna_ref = dna_ref_plus
        dna_alt = dna_ref[:snp_idx] + ref + dna_ref[snp_idx+1:]
        strand = "+ (alleles swapped: GRCh38 carries file-alt)"
    elif genomic_base_minus == ref:
        dna_ref = dna_ref_rc
        dna_alt = dna_ref[:snp_idx] + alt + dna_ref[snp_idx+1:]
        strand = "-"
    elif genomic_base_minus == alt:
        dna_ref = dna_ref_rc
        dna_alt = dna_ref[:snp_idx] + ref + dna_ref[snp_idx+1:]
        strand = "- (alleles swapped: GRCh38 carries file-alt on minus strand)"
    else:
        print(f"[WARN] {label}: allele mismatch. Plus={genomic_base}, "
              f"Minus={genomic_base_minus}, file ref={ref}, alt={alt}")
        dna_ref = dna_ref_plus
        dna_alt = dna_ref[:snp_idx] + alt + dna_ref[snp_idx+1:]
        strand = "? (unknown)"
    print(f"  Strand/orientation: {strand}")

    rna_ref = dna_to_rna(dna_ref)
    rna_alt = dna_to_rna(dna_alt)

    ref_hits = scan(rna_ref, snp_idx)
    alt_hits = scan(rna_alt, snp_idx)

    ref_set = {(r, m) for r, m, *_ in ref_hits}
    alt_set = {(r, m) for r, m, *_ in alt_hits}

    disrupted = [(r, m, ro) for r, m, ro, s, e in ref_hits if (r, m) not in alt_set]
    created   = [(r, m, ro) for r, m, ro, s, e in alt_hits if (r, m) not in ref_set]

    results[label] = {
        "gene": gene, "disease": disease,
        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
        "rna_ref_window": rna_ref[FLANK-10:FLANK+11],
        "rna_alt_window": rna_alt[FLANK-10:FLANK+11],
        "disrupted": disrupted,
        "created":   created,
        "all_ref_hits": ref_hits,
    }

# ── Report ────────────────────────────────────────────────────────────────────
for label, r in results.items():
    print(f"\n{'='*60}")
    print(f"{label}  ({r['chrom']}:{r['pos']}  {r['ref']}>{r['alt']})")
    print(f"  REF context: {r['rna_ref_window']}")
    print(f"  ALT context: {r['rna_alt_window']}")
    print(f"  All REF overlapping motifs: {len(r['all_ref_hits'])}")
    for rbp, iupac, role, s, e in r['all_ref_hits']:
        print(f"    {rbp:12s}  {iupac:15s}  [{role}]  window pos {s}-{e}")
    if r['disrupted']:
        print(f"  DISRUPTED by ALT allele ({len(r['disrupted'])}):")
        for rbp, motif, role in r['disrupted']:
            print(f"    *** {rbp:12s}  {motif:15s}  [{role}]")
    else:
        print("  No motifs disrupted.")
    if r['created']:
        print(f"  CREATED by ALT allele ({len(r['created'])}):")
        for rbp, motif, role in r['created']:
            print(f"    +++ {rbp:12s}  {motif:15s}  [{role}]")
