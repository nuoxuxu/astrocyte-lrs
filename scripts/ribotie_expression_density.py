import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import re


def load_ribotie_union(prefix):
    stim = pd.read_csv(f"{prefix}/ribotie_res_Stim.csv")
    unstim = pd.read_csv(f"{prefix}/ribotie_res_Unstim.csv")
    return set(stim["transcript_id"].unique()) | set(unstim["transcript_id"].unique())


def load_gtf_transcript_ids(gtf_path):
    tx_ids = set()
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] == "transcript":
                m = re.search(r'transcript_id "([^"]+)"', fields[8])
                if m:
                    tx_ids.add(m.group(1))
    return tx_ids


def prepare_df(ribotie_tx, expression_df):
    sample_cols = [c for c in expression_df.columns if c != "tname"]
    df = expression_df[["tname"]].copy()
    df["avg_count"] = expression_df[sample_cols].mean(axis=1)
    df["has_ribotie"] = df["tname"].isin(ribotie_tx)
    return df


def plot_hist(ax, df, title):
    detected = df.loc[df["has_ribotie"], "avg_count"]
    undetected = df.loc[~df["has_ribotie"], "avg_count"]

    # Use log scale for x; filter out zeros for log
    detected_log = np.log10(detected[detected > 0])
    undetected_log = np.log10(undetected[undetected > 0])

    hist_bins = np.linspace(-1, 4, 60)

    ax.hist(undetected_log, bins=hist_bins, alpha=0.5, color="#D65F5F",
            label=f"No RiboTIE (n={len(undetected):,})")
    ax.hist(detected_log, bins=hist_bins, alpha=0.5, color="#4C72B0",
            label=f"RiboTIE detected (n={len(detected):,})")

    ax.set_xlabel("Average read count (log10)")
    ax.set_ylabel("Number of transcripts")
    ax.set_title(title)
    ax.legend(fontsize=9)


# --- Load data ---
low_ribotie = load_ribotie_union("nextflow_results/ribotie/low_stringency")
low_orfanage_tx = load_gtf_transcript_ids("nextflow_results/orfanage/low_stringency/orfanage.gtf")
low_expr = pd.read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/low_stringency/final_expression.parquet"
)
low_expr = low_expr[low_expr["tname"].isin(low_orfanage_tx)]
low_df = prepare_df(low_ribotie, low_expr)

high_ribotie = load_ribotie_union("nextflow_results/ribotie/high_stringency")
high_orfanage_tx = load_gtf_transcript_ids("nextflow_results/orfanage/high_stringency/orfanage.gtf")
high_expr = pd.read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/high_stringency/final_expression.parquet"
)
high_expr = high_expr[high_expr["tname"].isin(high_orfanage_tx)]
high_df = prepare_df(high_ribotie, high_expr)

gencode_ribotie = load_ribotie_union("nextflow_results/ribotie/gencode")
low_clf = pd.read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/low_stringency/final_classification.parquet"
)
isoform_to_gencode = low_clf[["isoform", "associated_transcript"]].copy()
isoform_to_gencode = isoform_to_gencode[
    isoform_to_gencode["associated_transcript"].notna()
    & (isoform_to_gencode["associated_transcript"] != "novel")
]
sample_cols = [c for c in low_expr.columns if c != "tname"]
gencode_expr = low_expr.merge(
    isoform_to_gencode, left_on="tname", right_on="isoform", how="inner"
)
gencode_agg = gencode_expr.groupby("associated_transcript")[sample_cols].sum().reset_index()
gencode_agg.rename(columns={"associated_transcript": "tname"}, inplace=True)
gencode_df = prepare_df(gencode_ribotie, gencode_agg)

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

plot_hist(axes[0], low_df, "Low stringency")
plot_hist(axes[1], high_df, "High stringency")
plot_hist(axes[2], gencode_df, "GENCODE")

fig.suptitle("Distribution of average read counts\nby RiboTIE detection status",
             fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig("results/ribotie_expression_density.png", dpi=200, bbox_inches="tight")
plt.show()
