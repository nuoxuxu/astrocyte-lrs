import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import re

bins = [0, 1, 5, 10, 50, 100, np.inf]
labels = ["(0,1]", "(1,5]", "(5,10]", "(10,50]", "(50,100]", ">100"]


def compute_proportion(ribotie_tx, expression_df):
    sample_cols = [c for c in expression_df.columns if c != "tname"]
    expression_df = expression_df.copy()
    expression_df["avg_count"] = expression_df[sample_cols].mean(axis=1)

    df = expression_df[["tname", "avg_count"]].copy()
    df["has_ribotie"] = df["tname"].isin(ribotie_tx)
    df["count_bin"] = pd.cut(df["avg_count"], bins=bins, labels=labels, right=True)

    grouped = df.groupby("count_bin", observed=False).agg(
        total=("has_ribotie", "size"),
        with_ribotie=("has_ribotie", "sum"),
    )
    grouped["proportion"] = grouped["with_ribotie"] / grouped["total"]
    grouped = grouped[grouped["total"] > 0]
    return grouped


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


def plot_panel(ax, grouped, title):
    tick_labels = [str(idx) for idx in grouped.index]
    ax.bar(range(len(grouped)), grouped["proportion"], color="#4C72B0", edgecolor="white")
    for i, (idx, row) in enumerate(grouped.iterrows()):
        ax.text(i, row["proportion"] + 0.005, f"n={int(row['total'])}",
                ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(grouped)))
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Average read count across samples")
    ax.set_ylabel("Proportion with RiboTIE prediction")
    ax.set_title(title)
    ax.set_ylim(0, min(grouped["proportion"].max() * 1.15, 1.0))


# --- Low stringency ---
low_ribotie = load_ribotie_union("nextflow_results/ribotie/low_stringency")
low_orfanage_tx = load_gtf_transcript_ids("nextflow_results/orfanage/low_stringency/orfanage.gtf")
low_expr = pd.read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/low_stringency/final_expression.parquet"
)
low_expr = low_expr[low_expr["tname"].isin(low_orfanage_tx)]
low_grouped = compute_proportion(low_ribotie, low_expr)

# --- High stringency ---
high_ribotie = load_ribotie_union("nextflow_results/ribotie/high_stringency")
high_orfanage_tx = load_gtf_transcript_ids("nextflow_results/orfanage/high_stringency/orfanage.gtf")
high_expr = pd.read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/high_stringency/final_expression.parquet"
)
high_expr = high_expr[high_expr["tname"].isin(high_orfanage_tx)]
high_grouped = compute_proportion(high_ribotie, high_expr)

# --- GENCODE ---
gencode_ribotie = load_ribotie_union("nextflow_results/ribotie/gencode")

# Map long-read isoforms to GENCODE transcript IDs via classification
low_clf = pd.read_parquet(
    "nextflow_results/sqanti3/isoseq/sqanti3_filter/low_stringency/final_classification.parquet"
)
# Keep only isoforms that map to a real GENCODE transcript
isoform_to_gencode = low_clf[["isoform", "associated_transcript"]].copy()
isoform_to_gencode = isoform_to_gencode[
    isoform_to_gencode["associated_transcript"].notna()
    & (isoform_to_gencode["associated_transcript"] != "novel")
]

# Merge with expression to get per-isoform expression
gencode_expr = low_expr.merge(
    isoform_to_gencode, left_on="tname", right_on="isoform", how="inner"
)
sample_cols = [c for c in low_expr.columns if c != "tname"]

# Aggregate expression by GENCODE transcript (sum across isoforms mapping to same transcript)
gencode_agg = gencode_expr.groupby("associated_transcript")[sample_cols].sum().reset_index()
gencode_agg.rename(columns={"associated_transcript": "tname"}, inplace=True)

# Compute proportion using same binning logic
gencode_grouped = compute_proportion(gencode_ribotie, gencode_agg)

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

plot_panel(axes[0], low_grouped, "Low stringency")
plot_panel(axes[1], high_grouped, "High stringency")
plot_panel(axes[2], gencode_grouped, "GENCODE")

fig.suptitle("Proportion of transcripts with RiboTIE predictions\nby average expression level",
             fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig("results/ribotie_proportion_by_expression.png", dpi=200, bbox_inches="tight")
plt.show()

print("=== Low stringency ===")
print(low_grouped)
print("\n=== High stringency ===")
print(high_grouped)
print("\n=== GENCODE ===")
print(gencode_grouped)
