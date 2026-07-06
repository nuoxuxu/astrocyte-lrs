#!/scratch/nxu/astrocyte-lrs/env/bin/python3
"""Bar chart of in-frame ORF type distribution from orf_type_gencode.tsv.

Produces two separate figures — one for ORF_type_ORFanage and one for
ORF_type_RiboTIE. When --ribotie-gtf, --orfanage-gtf, and --gencode-gtf
are provided, each bar is overlaid with the fraction of ORFs in that category
whose CDS exon intervals exactly match a GENCODE v50 CDS (same chrom, strand,
and frozenset of (start, end) segments).
"""
import argparse
import csv
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_FRAME_CATEGORIES = [
    "annotated CDS",
    "N-terminal extension",
    "C-terminal extension",
    "N-terminal truncation",
    "C-terminal truncation",
]

_ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def _parse_attrs(attr_str):
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(attr_str)}


def load_cds_bins(gtf_path, key_attr):
    """Parse CDS features from a GTF, returning {key: (chrom, strand, frozenset{(start,end)})}."""
    segs = defaultdict(set)
    meta = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = _parse_attrs(fields[8])
            key = attrs.get(key_attr)
            if not key:
                continue
            meta[key] = (fields[0], fields[6])
            segs[key].add((int(fields[3]), int(fields[4])))
    return {k: (*meta[k], frozenset(segs[k])) for k in segs}


def load_gencode_cds_bin_set(gtf_path):
    """Return a set of (chrom, strand, frozenset{(start,end)}) for every GENCODE CDS transcript."""
    segs = defaultdict(set)
    meta = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = _parse_attrs(fields[8])
            tx = attrs.get('transcript_id', '').split('.')[0]
            if not tx:
                continue
            meta[tx] = (fields[0], fields[6])
            segs[tx].add((int(fields[3]), int(fields[4])))
    return {(*meta[k], frozenset(segs[k])) for k in segs}


def count_gencode_matches(keys_by_cat, cds_bins, gencode_bin_set):
    """For each category, count how many keys have an exact GENCODE CDS bin match."""
    match_counts = Counter()
    for cat, keys in keys_by_cat.items():
        for key in keys:
            entry = cds_bins.get(key)
            if entry and entry in gencode_bin_set:
                match_counts[cat] += 1
    return match_counts


def plot_one(counts, total, match_counts, column_label, version_name, out_path):
    cats = IN_FRAME_CATEGORIES
    vals = [counts.get(c, 0) for c in cats]
    match_vals = [match_counts.get(c, 0) if match_counts is not None else 0 for c in cats]
    max_val = max(vals) if any(vals) else 1
    bar_h = 0.55

    fig, ax = plt.subplots(figsize=(7.5, 4))

    # Full bar: light shade
    ax.barh(cats, vals, height=bar_h, color="#1565C0", alpha=0.35, edgecolor="none")
    # Overlay: exact GENCODE CDS match, darker shade
    if match_counts is not None:
        ax.barh(cats, match_vals, height=bar_h, color="#1565C0", alpha=0.90,
                edgecolor="none", label="exact GENCODE CDS match")

    offset = max_val * 0.012
    for i, (val, mval) in enumerate(zip(vals, match_vals)):
        pct = 100 * val / total if total else 0
        # Count + % of in-frame total, right of full bar
        ax.text(val + offset, i, f"{val:,}  ({pct:.1f}%)",
                va="center", ha="left", fontsize=9)
        # % matching GENCODE, inside the overlay bar in white
        if match_counts is not None and val > 0:
            mpct = 100 * mval / val
            if mval > max_val * 0.04:
                ax.text(mval - offset, i, f"{mpct:.0f}%",
                        va="center", ha="right", fontsize=8,
                        color="white", fontweight="bold")

    suffix = f" — {version_name}" if version_name else ""
    ax.set_title(f"{column_label}{suffix}", fontsize=11)
    ax.set_xlabel("Number of ORFs", fontsize=10)
    if match_counts is not None:
        ax.legend(fontsize=8, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, max_val * 1.40)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("orf_type_tsv", help="orf_type_gencode.tsv")
    parser.add_argument("--name", default="", help="Version label (e.g. no_minlen)")
    parser.add_argument("--prefix", default="orf_type_bars", help="Output filename prefix")
    parser.add_argument("--ribotie-gtf",  help="RiboTIE merged GTF (CDS keyed by ORF_id)")
    parser.add_argument("--orfanage-gtf", help="Plain orfanage GTF (CDS keyed by transcript_id)")
    parser.add_argument("--gencode-gtf",  help="GENCODE annotation GTF for exact CDS bin matching")
    args = parser.parse_args()

    orfanage_counts   = Counter()
    ribotie_counts    = Counter()
    orfanage_keys_by_cat = defaultdict(list)   # cat -> [transcript_id, ...]
    ribotie_keys_by_cat  = defaultdict(list)   # cat -> [ORF_id, ...]

    with open(args.orf_type_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            oa, rt = row["ORF_type_ORFanage"], row["ORF_type_RiboTIE"]
            tx_id, orf_id = row["transcript_id"], row["ORF_id"]
            if oa in IN_FRAME_CATEGORIES:
                orfanage_counts[oa] += 1
                orfanage_keys_by_cat[oa].append(tx_id)
            if rt in IN_FRAME_CATEGORIES:
                ribotie_counts[rt] += 1
                ribotie_keys_by_cat[rt].append(orf_id)

    oa_match = None
    rt_match = None

    if args.ribotie_gtf and args.orfanage_gtf and args.gencode_gtf:
        print("Loading GENCODE CDS bins...")
        gencode_bins = load_gencode_cds_bin_set(args.gencode_gtf)
        print(f"  {len(gencode_bins):,} GENCODE CDS bins loaded")

        oa_cds = load_cds_bins(args.orfanage_gtf, "transcript_id")
        oa_match = count_gencode_matches(orfanage_keys_by_cat, oa_cds, gencode_bins)

        rt_cds = load_cds_bins(args.ribotie_gtf, "ORF_id")
        rt_match = count_gencode_matches(ribotie_keys_by_cat, rt_cds, gencode_bins)

    plot_one(orfanage_counts, sum(orfanage_counts.values()), oa_match,
             "ORF_type_ORFanage", args.name, f"{args.prefix}_orfanage.png")
    plot_one(ribotie_counts, sum(ribotie_counts.values()), rt_match,
             "ORF_type_RiboTIE", args.name, f"{args.prefix}_ribotie.png")


if __name__ == "__main__":
    main()
