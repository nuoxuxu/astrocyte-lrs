#!/bin/env python
import polars as pl
import argparse

# IsoformSwitchAnalyzeR columns (condition_1 = Unstim, condition_2 = Stim; dIF = IF2 - IF1)
ISOSWITCH_COLS = [
    "gene_log2_fold_change", "gene_q_value",
    "iso_log2_fold_change", "iso_q_value",
    "IF1", "IF2", "dIF",
    "isoform_switch_q_value", "gene_switch_q_value",
    "switchConsequencesGene",
]

# limma Stim (+) vs Unstim (-) protein differential expression columns
PROTEIN_DE_COLS = [
    "log2FC", "AveExpr", "P.Value", "adj.P.Val",
    "n_valid_Unstim", "n_valid_Stim", "n_imputed", "direction",
]


def main():
    parser = argparse.ArgumentParser(
        description='Add differential splicing (IsoformSwitchAnalyzeR) and differential protein '
                    'expression (limma) columns to the RiboTIE ORF table, joined on ORF_id.')
    parser.add_argument('ribotie_csv', help='RiboTIE ORF CSV (ribotie_res_merged_fixed_with_lncRNA.csv)')
    parser.add_argument('isoform_features', help='IsoformSwitchAnalyzeR isoformFeatures.csv keyed by ORF_id')
    parser.add_argument('protein_de', help='protein_DE_stim_vs_unstim.tsv')
    parser.add_argument('-o', '--output', default='ribotie_summary.csv', help='Output CSV file')
    args = parser.parse_args()

    ribotie = pl.read_csv(args.ribotie_csv)

    isoswitch = pl.read_csv(args.isoform_features, null_values="NA", infer_schema_length=None)\
        .select(
            pl.col("isoform_id").alias("ORF_id"),
            *[pl.col(c).alias(f"isoswitch_{c}") for c in ISOSWITCH_COLS]
        )\
        .rename({"isoswitch_IF1": "isoswitch_IF_Unstim", "isoswitch_IF2": "isoswitch_IF_Stim"})

    protein_de = pl.read_csv(args.protein_de, separator="\t", infer_schema_length=None)\
        .select(
            pl.col("Protein.Group").alias("ORF_id"),
            *[pl.col(c).alias(f"protein_{c.replace('.', '_')}") for c in PROTEIN_DE_COLS]
        )

    ribotie\
        .join(isoswitch, on="ORF_id", how="left")\
        .join(protein_de, on="ORF_id", how="left")\
        .write_csv(args.output)


if __name__ == "__main__":
    main()
