import polars as pl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

df = pl.read_csv(
    "nextflow_results/quality/minlen_quality_metrics_riboseq.tsv",
    separator="\t",
    null_values=["NA", "UnableToDownloadAlignment"],
)

novel = df.filter(pl.col("is_novel") == 1)["ORBLv"].drop_nulls().to_numpy()
known = df.filter(pl.col("is_novel") == 0)["ORBLv"].drop_nulls().to_numpy()

fig, ax = plt.subplots(figsize=(7, 4))

bins = 50
ax.hist(known, bins=bins, alpha=0.6, label=f"Known (n={len(known):,})", color="#4878CF", density=True)
ax.hist(novel, bins=bins, alpha=0.6, label=f"Novel (n={len(novel):,})", color="#D65F5F", density=True)

ax.set_xlabel("ORBLv score")
ax.set_ylabel("Density")
ax.set_title("ORBLv score distribution by novelty (minlen)")
ax.legend()
ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())

plt.tight_layout()
plt.savefig("scripts/orblv_by_is_novel.png", dpi=150)
print("Saved scripts/orblv_by_is_novel.png")
