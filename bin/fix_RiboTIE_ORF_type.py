#!/usr/bin/env python
import os
import sys
import argparse
import polars as pl
import h5py
from tqdm import tqdm
import numpy as np
from transcript_transformer.util_functions import prtime
from transcript_transformer.processing import (
    parse_CDS_overlap, 
    save_output_table, 
    filter_CDS_variants
)

try:
    from transcript_transformer.util_functions import load_config
except ImportError:
    # Installed transcript_transformer (e.g. astrocytes/pytorch venv) predates load_config;
    # copied from transcript_transformer/util_functions.py
    def load_config(config_path, tool="tis_transformer"):
        """Load configuration from a YAML/JSON file for interactive use.

        Args:
            config_path (str): Path to the configuration file (YAML or JSON)
            tool (str): Tool name, either "tis_transformer" or "ribotie". Defaults to "tis_transformer"

        Returns:
            argparse.Namespace: Arguments object with config parameters
        """
        from transcript_transformer.argparser import Parser
        from importlib.resources import files
        from typing import cast

        parser = Parser(stage="train", tool=tool)
        run_parser = parser.add_run_args()

        if tool == "ribotie":
            run_parser.add_argument(
                "--pretrain",
                action="store_true",
                help="pretrain model using all available samples.",
            )

        data_parser = parser.add_data_args()

        if tool == "tis_transformer":
            data_parser.add_argument(
                "--model",
                choices=["human", "mouse"],
                default=None,
                type=str,
                help="Use a pre-trained model for predictions.",
            )
            fa_parser = parser.add_argument_group(
                "fasta",
                "Use a fasta file with sequences to predict.",
            )
            fa_parser.add_argument(
                "--fasta",
                type=str,
                default=None,
                help="Path to a fasta file with sequences to predict.",
            )
            fa_parser.add_argument(
                "--fold",
                type=int,
                default=0,
                help="Determines the fold number to use for prediction.",
            )

        pr_parser = parser.add_processing_args()
        if tool == "ribotie":
            pr_parser.add_argument(
                "--no_correction",
                action="store_true",
                help="Don't correct to nearest in-frame ATG.",
            )
            pr_parser.add_argument(
                "--distance",
                type=int,
                default=9,
                help="Number of codons to search up- and downstream for an ATG.",
            )

        parser.add_comp_args()
        parser.add_training_args()
        parser.add_train_loading_args()
        parser.add_evaluation_args()

        # Load default config if available
        try:
            if tool == "tis_transformer":
                default_config = files("transcript_transformer.configs").joinpath("defaults.tt.yml")
            else:
                default_config = files("transcript_transformer.configs").joinpath("defaults.rt.yml")
            default_config = os.fspath(cast(os.PathLike, default_config))
            args = parser.parse_arguments([config_path], configs=[default_config, config_path])
        except Exception:
            # If default config not found, just load the provided config
            args = parser.parse_arguments([config_path], configs=[config_path])

        return args

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                ).drop("attributes")

def get_df_CDS(h5_path):
    # Taken from construct_output_table
    # load in all CDS properties in h5 db
    f = h5py.File(h5_path, "r")
    h5_cols = [
        "transcript_id",
        "seqname",
        "strand",
        "CDS_coords",
        "canonical_TIS_coord",
        "canonical_LTS_coord",
    ]
    mask = pl.Series(list(f[f"transcript/canonical_TIS_idx"])) != -1
    df_CDS = (
        pl.DataFrame(
            {f"{h}": np.array(f[f"transcript/{h}"])[mask.arg_true()] for h in h5_cols}
        )
        .with_columns(
            pl.col("CDS_coords").map_elements(list, pl.List(pl.Int64)),
            pl.col(pl.Binary).cast(pl.String),
        )
        .with_columns(
            CDS_exon_start=pl.col("CDS_coords").list.gather_every(2, 0),
            CDS_exon_end=pl.col("CDS_coords").list.gather_every(2, 1),
            CDS_start_range=(
                pl.when(pl.col("strand") == "+")
                .then(pl.col("CDS_coords").list.get(0))
                .otherwise(pl.col("CDS_coords").list.get(-2))
            ),
            CDS_end_range=(
                pl.when(pl.col("strand") == "+")
                .then(pl.col("CDS_coords").list.get(-1))
                .otherwise(pl.col("CDS_coords").list.get(1))
            ),
        )
        .drop("CDS_coords")
    )
    # close h5 db handle
    f.file.close()
    return df_CDS

out_headers = ['ribotie_score', 'ribotie_rank', 'seqname', 'ORF_id', 'ORF_len', 'transcript_id', 'transcript_len', 'start_codon', 'stop_codon', 'strand', 'ORF_type', 
               'TIS_pos', 'TTS_pos', 'has_CDS_clones', 'has_CDS_TIS', 'has_CDS_TTS', 'shared_in_frame_CDS_frac', 'dist_from_canonical_TIS', 'frame_wrt_canonical_TIS', 
               'TTS_on_transcript', 'TIS_coord', 'TIS_exon', 'TTS_coord', 'TTS_exon', 'LTS_coord', 'LTS_exon', 'gene_id', 'canonical_TIS_coord', 'canonical_TIS_pos', 
               'canonical_LTS_coord', 'canonical_LTS_pos', 'canonical_TTS_coord', 'canonical_TTS_pos', 'has_annotated_start_codon', 'has_annotated_stop_codon', 
               'protein_seq', 'correction', 'reads_in_transcript', 'reads_in_ORF', 'reads_in_frame_frac', 'reads_5UTR', 'reads_3UTR', 'reads_coverage_frac', 
               'reads_entropy', 'reads_skew', 'gene_name']

# Load RiboTIE output


def build_df(gtf_path, ribotie_cpm1_3sample_path):
    """
    Build df for parse_CDS_overlap and filter_CDS_variants
    
    :param gtf: Description
    :param ribotie_cpm1_3sample: Description
    """
    # Step 1: Extract ORF CDS exons from GTF and sort within each ORF
    # For + strand: sort start ascending; for - strand: sort start descending
    ribotie_cpm1_3sample = pl.read_csv(ribotie_cpm1_3sample_path)
    gtf = read_gtf(gtf_path)

    df_cds = gtf\
            .filter(pl.col("feature") == "CDS")\
            .select("transcript_id", "start", "end", "strand")\
            .with_columns(
                sort_key=pl.when(pl.col("strand") == "+")
                .then(pl.col("start"))
                .otherwise(-pl.col("start"))
            )\
            .sort(["transcript_id", "sort_key"])\
            .drop("sort_key")

    # Step 2: Aggregate CDS regions per ORF
    df_cds = (
        df_cds
        .group_by("transcript_id", maintain_order=True)
        .agg(
            pl.col("start").alias("ORF_exon_start"),
            pl.col("end").alias("ORF_exon_end"),
            pl.col("strand").first()
        )
    )

    # Compute exon lengths: (end - start) + 1 for each exon pair
    df_cds = df_cds.with_columns(
        ORF_exon_len=(
            pl.col("ORF_exon_end") - pl.col("ORF_exon_start") + 1
        )
    ).rename({"transcript_id": "ORF_id"})

    # Step 3: Prepare RiboTIE CSV: drop problematic columns that will be overwritten by parse_CDS_overlap
    cols_to_drop = ["", "has_CDS_clones", "has_CDS_TIS", "has_CDS_TTS", "shared_in_frame_CDS_frac"]
    cols_to_drop = [c for c in cols_to_drop if c in ribotie_cpm1_3sample.columns]
    ribotie_prepared = ribotie_cpm1_3sample.drop(cols_to_drop)

    # Step 4: Join GTF CDS exons with RiboTIE data
    df = ribotie_prepared.join(
        df_cds,
        left_on="ORF_id",
        right_on="ORF_id",
        how="left"
    )
    return df

def save_output_table(df, out_prefix, label, prefix, out_headers):
    df = (
        df.with_columns(
            (pl.col(f"{prefix}score").rank(method="ordinal", descending=True)).alias(
                f"{prefix}rank"
            )
        )
        .select(out_headers)
        .sort(f"{prefix}rank")
    )
    df.write_csv(f"{out_prefix}{label}.csv", float_precision=4)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge RiboTIE ORF predictions with CDS annotations from RiboTIE h5 databases, "
        "re-evaluate CDS overlap/variants and write redundant, filtered and novel ORF tables."
    )
    parser.add_argument("--config", required=True, help="RiboTIE YAML config (provides start_codons, min_ORF_len, include_invalid_TTS)")
    parser.add_argument("--gtf", required=True, help="GTF of predicted ORFs (transcript_id == ORF_id)")
    parser.add_argument("--ribotie_csv", required=True, help="RiboTIE ORF predictions CSV")
    parser.add_argument("--h5", required=True, nargs="+", help="RiboTIE h5 database(s) whose annotated CDSs are merged")
    parser.add_argument("--out_prefix", required=True, help="Output prefix; writes {prefix}.redundant.csv, {prefix}.csv and {prefix}.novel.csv")
    return parser.parse_args()

def main(cli):
    # Parser.parse_arguments re-parses sys.argv; hide this script's CLI flags from it
    argv, sys.argv = sys.argv, sys.argv[:1]
    args = load_config(cli.config, tool="ribotie")
    sys.argv = argv
    df_ORFanage = build_df(
        gtf_path=cli.gtf,
        ribotie_cpm1_3sample_path=cli.ribotie_csv
    )

    # Maybe we can just combine GENCODE and ORFanage annotation as input to RiboTIE
    df_CDS_merged = pl.concat([get_df_CDS(h5_path=h5_path) for h5_path in cli.h5])

    # To evaluate CDS variants, group df and df_CDS by seqname (to prevent OOM)
    df_grps = []
    total = df_ORFanage["seqname"].unique().len()
    for seqname, df_grp in tqdm(df_ORFanage.group_by("seqname"), total=total, desc="seqname"):
        df_CDS_grp = df_CDS_merged.filter(pl.col("seqname") == seqname[0])
        df_grp = parse_CDS_overlap(df_grp, df_CDS_grp)
        df_grps.append(df_grp)

    df = pl.concat(df_grps)
    df = df.with_columns(pl.col("shared_in_frame_CDS_frac").truediv(pl.col("ORF_len")))

    # --- Filter CDS variants and custom filters ---
    # Custom filters
    n_total = len(df)
    n_before = n_total
    if not args.include_invalid_TTS:
        df = df.filter(pl.col("TTS_on_transcript"))
        n_after = len(df)
        n_rem = n_before - n_after
        prtime(
            f"{n_rem} of {n_total} ({n_rem/n_total:.1%}) ORFs removed due to invalid TTS",
            "\t",
        )
        n_before = n_after

    df = df.filter(pl.col("start_codon").str.contains(args.start_codons))
    n_after = len(df)
    n_rem = n_before - n_after
    prtime(
    f"{n_rem} of {n_total} ({n_rem/n_total:.1%}) ORFs removed due to start codon filter",
    "\t",
    )
    n_before = n_after

    df = df.filter(pl.col("ORF_len") >= args.min_ORF_len)
    n_after = len(df)
    n_rem = n_before - n_after
    prtime(
    f"{n_rem} of {n_total} ({n_rem/n_total:.1%}) ORFs removed due to minimum ORF length filter",
    "\t",
    )

    # CDS variant filtering
    if len(df) > 0:
        df_filt = filter_CDS_variants(df)
    else:
        df_filt = df
    df_novel = df_filt.filter(pl.col("ORF_type") != "annotated CDS")

    # --- Save to csv ---
    prefix = "ribotie_"
    out_headers.remove("has_annotated_start_codon")
    out_headers.remove("has_annotated_stop_codon")
    for df_, label in zip([df, df_filt, df_novel], [".redundant", "", ".novel"]):
        save_output_table(df_, cli.out_prefix, label, prefix, out_headers)


if __name__ == "__main__":
    main(parse_args())
