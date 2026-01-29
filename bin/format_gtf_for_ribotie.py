#!/bin/env python
import polars as pl
from src.utils import read_gtf
import polars as pl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_gtf", action="store", dest="input_gtf")
parser.add_argument("--final_classification", action="store", dest="final_classification")
parser.add_argument("--annotation_gtf", action="store", dest="annotation_gtf")
parser.add_argument("--output_gtf", action="store", dest="output_gtf")
params = parser.parse_args()

def add_gene_name(df, classification, annotation_gtf):
    import polars as pl
    id_to_gene_name = read_gtf(annotation_gtf, attributes=["gene_id", "gene_name"])\
        .filter(pl.col("feature")=="gene")\
        .select("gene_id", "gene_name")\
        .unique()
    classification = pl.read_parquet(classification)\
        .join(
            id_to_gene_name,
            left_on="associated_gene",
            right_on="gene_id",
            how="left"
        )
    return df.join(
        classification.select(["isoform", "gene_name"]),
        left_on="transcript_id", 
        right_on="isoform"
    )

def write_gtf(df: pl.DataFrame, output_path: str):
    std_cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame"]
    
    # 1. Identify attribute columns
    attr_cols = [c for c in df.columns if c not in std_cols]
    
    # 2. Define which attributes should NOT be quoted
    # Add any other numeric attributes here if needed
    unquoted_attrs = {"exon_number"} 

    attr_exprs = []
    for col_name in attr_cols:
        # Determine the format string based on whether the column is in our unquoted list
        if col_name in unquoted_attrs:
            fmt = '{} {};'   # Result: exon_number 1;
        else:
            fmt = '{} "{}";' # Result: gene_id "ENSG...";
            
        expr = (
            pl.when(pl.col(col_name).is_not_null())
            .then(pl.format(fmt, pl.lit(col_name), pl.col(col_name)))
            .otherwise(None)
        )
        attr_exprs.append(expr)

    # 3. Construct the export DataFrame
    export_df = df.with_columns([
        pl.col("score").cast(pl.String).fill_null("."),
        pl.col("frame").cast(pl.String).fill_null("."),
        
        # Concatenate attributes
        pl.concat_str(attr_exprs, separator=" ", ignore_nulls=True).alias("attributes")
    ]).select(std_cols + ["attributes"])

    # 4. Write to disk
    export_df.write_csv(
        output_path, 
        separator="\t", 
        include_header=False, 
        quote_style="never"
    )
    print(f"Written to {output_path}")

attribute_list = ["gene_id", "transcript_id"]
orfanage_gtf = read_gtf(params.input_gtf, attributes=attribute_list)

orfanage_exons = orfanage_gtf.filter(pl.col("feature") == "exon")
orfanage_CDS = orfanage_gtf.filter(pl.col("feature") == "CDS")

exons_numbered_df = orfanage_exons\
    .with_columns(
        exon_number = pl.when(pl.col("strand") == "-")
        .then(
            # For negative strand, the highest start position is Exon 1
            pl.col("start").rank(method="ordinal", descending=True).over("transcript_id")
        )
        .otherwise(
            # For positive strand, the lowest start position is Exon 1
            pl.col("start").rank(method="ordinal", descending=False).over("transcript_id")
        )
    )

orfanage_CDS = orfanage_CDS\
    .join_where(
        exons_numbered_df.select(["transcript_id", "start", "end", "exon_number"]),
        pl.col("start") >= pl.col("start_right"),
        pl.col("end") <= pl.col("end_right"),
        pl.col("transcript_id") == pl.col("transcript_id_right")
    ).drop(["transcript_id_right", "start_right", "end_right"])

numbered_gtf = pl.concat([exons_numbered_df, orfanage_CDS], how="vertical").sort(["transcript_id", "start"])
other_feature = orfanage_gtf.filter(pl.col("feature").is_in(["exon", "CDS"]).not_())\
    .with_columns(exon_number = pl.lit(None))\
    .with_columns(pl.col("exon_number").cast(pl.UInt32))

out_gtf = pl.concat([other_feature, numbered_gtf], how="vertical")\
    .pipe(add_gene_name, params.final_classification, params.annotation_gtf)\
    .drop("attributes")\
    .sort(["gene_id", "transcript_id", "start", "feature"])

write_gtf(out_gtf, params.output_gtf)