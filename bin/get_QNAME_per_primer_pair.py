import polars as pl
from pathlib import Path

# Get the list of unique primer pairs
primer_pair_list = list(set([str(path).split(".")[1] for path in Path("nextflow_results/prepare/refine").glob("*.flnc.report.csv")]))

# Get the list of flnc.report.csv files for each primer pair
def get_files_containing_primer_pair(primer_pair):
    return [path for path in list(Path("nextflow_results/prepare/refine").glob("*.flnc.report.csv")) if "IsoSeqX_bc08_5p--IsoSeqX_3p" in str(path)]

# Write all QNAMEs corresponding to each primer pair to a list
for primer_pair in primer_pair_list:
    pl.scan_csv(get_files_containing_primer_pair(primer_pair)).select("id").collect().write_csv(f"proc/{primer_pair}_QNAME_list.csv", include_header=False)

# Get the read_stat file from the merged collaped transcriptome
merged_collapsed_read_stat = pl.read_csv("nextflow_results/discover/isoseq/merged.collapsed.read_stat.txt", separator="\t")

# Testing one primer pair only
IsoSeqX_bc07_5p_list = pl.scan_csv(get_files_containing_primer_pair("IsoSeqX_bc07_5p--IsoSeqX_3p")).collect()["id"]

merged_collapsed_read_stat\
    .with_columns(
        pl.col("id").is_in(IsoSeqX_bc07_5p_list.to_list()).alias("IsoSeqX_bc07_5p--IsoSeqX_3p")
    )

# Gives segementfault error, probably out of memory