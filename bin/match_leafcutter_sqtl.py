#!/usr/bin/env python3
"""
Match significant leafcutter intron clusters with BigBrain cis-sQTLs.

Steps:
  1. Filter leafcutter_ds_cluster_significance.txt: status=="Success" and p.adjust<0.05
  2. Extract cluster names (strip chr prefix), match against leafcutter.clu2gene.txt col 1
  3. Subtract 1 from end coordinate; build feature string chr_start_end
  4. Join with BigBrain_cis_sQTL_ALL_top_assoc.tsv on feature; keep meta-analysis columns
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


def main():
    parser = argparse.ArgumentParser(
        description="Match significant leafcutter intron clusters with BigBrain cis-sQTLs."
    )
    parser.add_argument("sig_file",   help="leafcutter_ds_cluster_significance.txt")
    parser.add_argument("clu2gene",   help="leafcutter.clu2gene.txt")
    parser.add_argument("sqtl_file",  help="BigBrain cis-sQTL TSV file")
    parser.add_argument("-o", "--output", default="leafcutter_sqtl_matches.tsv",
                        help="Output TSV file (default: leafcutter_sqtl_matches.tsv)")
    args = parser.parse_args()

    sig = pl.read_csv(args.sig_file, separator="\t", null_values="NA")
    sig_filtered = sig.filter(
        (pl.col("status") == "Success") & (pl.col("p.adjust") < 0.05)
    )
    cluster_names = (
        sig_filtered["cluster"]
        .str.split(":")
        .list.get(1)
        .to_list()
    )
    print(f"Significant clusters (p.adjust < 0.05): {len(cluster_names)}")

    clu2gene = pl.read_csv(
        args.clu2gene,
        separator=" ",
        has_header=False,
        new_columns=["clu_name", "chr", "start", "end", "gene_id", "m1", "m2", "m3", "m4"],
    )
    matched = clu2gene.filter(pl.col("clu_name").is_in(cluster_names))
    print(f"clu2gene rows matched: {len(matched)}")

    matched = matched.with_columns((pl.col("end") - 1).alias("end"))
    matched = matched.with_columns(
        (
            pl.col("chr") + "_"
            + pl.col("start").cast(pl.String) + "_"
            + pl.col("end").cast(pl.String)
        ).alias("feature")
    )
    features = set(matched["feature"].to_list())

    sqtl = pl.read_csv(args.sqtl_file, separator="\t", null_values="NA")
    sqtl_filtered = sqtl.filter(pl.col("feature").is_in(features))
    keep_cols = [c for c in META_COLS if c in sqtl_filtered.columns]
    result = sqtl_filtered.select(keep_cols)
    print(f"sQTL rows matched: {len(result)}")

    result.write_csv(args.output, separator="\t")
    print(f"Output written to: {args.output}")


if __name__ == "__main__":
    main()
