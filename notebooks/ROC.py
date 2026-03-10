import marimo

__generated_with = "0.20.4"
app = marimo.App()


@app.cell
def _():
    import marimo as mo
    from gtfparse import read_gtf
    import polars as pl
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).parent.parent))
    from bin.src.trim_ribotish_stop_codon import trim_stop_codon
    from bin.src.label_orfs_by_sequence import label_noncanonical_orfs_by_sequence, extract_gtf_cds_proteins
    from bin.src.roc_data import label_query_transcripts, pct_annotated, build_roc_dfs

    return (
        build_roc_dfs,
        extract_gtf_cds_proteins,
        label_noncanonical_orfs_by_sequence,
        label_query_transcripts,
        mo,
        np,
        pct_annotated,
        pl,
        plt,
        read_gtf,
        roc_auc_score,
        roc_curve,
        trim_stop_codon,
    )


@app.cell
def _(read_gtf):
    ribotish_gtf = read_gtf("nextflow_results/ribotish/mid_stringency/ORFannotate_annotated.gtf")
    ribotie_gtf = read_gtf("nextflow_results/ribotie/mid_stringency/merged/ribotie_res_merged.redundant.gtf")
    gencode_gtf = read_gtf("data/gencode.v47.annotation.gtf")
    return gencode_gtf, ribotie_gtf, ribotish_gtf


@app.cell
def _(pl):
    ORFannotate_summary = pl.read_csv(
        "nextflow_results/ribotish/mid_stringency/ORFannotate_summary.tsv", separator="\t"
    )
    return (ORFannotate_summary,)


@app.cell
def _(pl, ribotish_gtf, trim_stop_codon):
    n_before = ribotish_gtf.filter(pl.col("feature") == "CDS").height
    trimmed = trim_stop_codon(ribotish_gtf)
    n_after = trimmed.filter(pl.col("feature") == "CDS").height
    dropped = n_before - n_after
    if dropped:
        print(f"Dropped {dropped} CDS rows that were entirely stop codon")
    print(f"CDS rows: {n_before} → {n_after}")
    return (trimmed,)


@app.cell
def _(pl):
    transcode_phase1 = pl.read_excel(
        "data/41587_2022_1369_MOESM2_ESM.xlsx", sheet_id=3, infer_schema_length=None
    )
    transcode_phase2 = pl.read_excel(
        "data/media-1.xlsx", sheet_id=3, infer_schema_length=None
    )
    return transcode_phase1, transcode_phase2


@app.cell
def _(mo):
    mo.md("""
    # ROC vs GENCODE
    """)
    return


@app.cell
def _(gencode_gtf, label_query_transcripts, ribotie_gtf, trimmed):
    ribotish_gtf_annotation = label_query_transcripts(trimmed, gencode_gtf)
    ribotie_gtf_annotation = label_query_transcripts(ribotie_gtf, gencode_gtf)
    return ribotie_gtf_annotation, ribotish_gtf_annotation


@app.cell
def _(pl, ribotie_gtf, ribotie_gtf_annotation):
    ribotie_scores = (
        ribotie_gtf.filter(pl.col("feature") == "CDS").unique("transcript_id")
        .join(ribotie_gtf_annotation, on="transcript_id", how="left", coalesce=True)
        .select(["ribotie_score", "label"])
        .with_columns(pl.col("ribotie_score").cast(pl.Float64))
        .with_columns(pl.col("label").replace({"annotated": 1, "novel": 0}).cast(pl.Int64))
    )
    return (ribotie_scores,)


@app.cell
def _(ORFannotate_summary, np, pl, ribotish_gtf_annotation):
    ribotish_scores = (
        ORFannotate_summary
        .join(ribotish_gtf_annotation, on="transcript_id", how="left", coalesce=True)
        .select(["frame_qvalue", "label"])
        .drop_nulls()
        .with_columns(
            (-np.log10(pl.col("frame_qvalue").cast(pl.Float64))).alias("frame_qvalue"),
            pl.col("label").replace({"annotated": 1, "novel": 0}).cast(pl.Int64),
        )
    )
    return (ribotish_scores,)


@app.cell
def _(plt, ribotie_scores, ribotish_scores, roc_auc_score, roc_curve):
    rbt_fpr, rbt_tpr, _ = roc_curve(ribotie_scores["label"], ribotie_scores["ribotie_score"], pos_label=1)
    rbt_auc = roc_auc_score(ribotie_scores["label"], ribotie_scores["ribotie_score"])
    rbtish_fpr, rbtish_tpr, _ = roc_curve(ribotish_scores["label"], ribotish_scores["frame_qvalue"], pos_label=1)
    rbtish_auc = roc_auc_score(ribotish_scores["label"], ribotish_scores["frame_qvalue"])

    plt.figure()
    plt.plot(rbt_fpr, rbt_tpr, label=f"ORFanage + RiboTIE (AUC = {rbt_auc:.2f})")
    plt.plot(rbtish_fpr, rbtish_tpr, label=f"Ribo-TISH (AUC = {rbtish_auc:.2f})")
    plt.plot([0, 1], [0, 1], "k--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC vs GENCODE")
    plt.legend()
    plt.show()
    return


@app.cell
def _(mo):
    mo.md("""
    # TransCode Proteomics Validation
    """)
    return


@app.cell
def _(transcode_phase1, transcode_phase2):
    GENOME_FASTA = "/Users/xunuo/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa"
    phase1_seqs = {str(s).rstrip("*") for s in transcode_phase1["orf_sequence"] if s is not None}
    phase2_seqs = {str(s).rstrip("*") for s in transcode_phase2["sequence_aa_MS"] if s is not None}
    return GENOME_FASTA, phase1_seqs, phase2_seqs


@app.cell
def _(
    GENOME_FASTA,
    label_noncanonical_orfs_by_sequence,
    ribotie_gtf,
    transcode_phase1,
    trimmed,
):
    ribotish_phase_1 = label_noncanonical_orfs_by_sequence(
        trimmed, transcode_phase1, GENOME_FASTA
    ).select("label")
    ribotie_phase_1 = label_noncanonical_orfs_by_sequence(
        ribotie_gtf, transcode_phase1, GENOME_FASTA
    ).select("label")
    return ribotie_phase_1, ribotish_phase_1


@app.cell
def _(
    GENOME_FASTA,
    label_noncanonical_orfs_by_sequence,
    ribotie_gtf,
    transcode_phase2,
    trimmed,
):
    ribotish_phase_2 = label_noncanonical_orfs_by_sequence(
        trimmed, transcode_phase2, GENOME_FASTA, sequence_col="sequence_aa_MS"
    ).select("label")
    ribotie_phase_2 = label_noncanonical_orfs_by_sequence(
        ribotie_gtf, transcode_phase2, GENOME_FASTA, sequence_col="sequence_aa_MS"
    ).select("label")
    return ribotie_phase_2, ribotish_phase_2


@app.cell
def _(
    pct_annotated,
    plt,
    ribotie_phase_1,
    ribotie_phase_2,
    ribotish_phase_1,
    ribotish_phase_2,
):
    tools = ["Ribo-TISH", "RiboTIE"]
    x = range(len(tools))
    width = 0.5

    fig, axes = plt.subplots(1, 2, figsize=(8, 5), sharey=True)

    for ax, rbtish, rbt, phase in [
        (axes[0], ribotish_phase_1, ribotie_phase_1, "Phase 1"),
        (axes[1], ribotish_phase_2, ribotie_phase_2, "Phase 2"),
    ]:
        vals = [pct_annotated(rbtish), pct_annotated(rbt)]
        bars = ax.bar(list(x), vals, width)
        ax.set_xticks(list(x))
        ax.set_xticklabels(tools)
        ax.set_title(phase)
        ax.set_ylabel("% ORFs annotated")
        ax.set_ylim(0, 100)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, val + 1, f"{val:.1f}%", ha="center", va="bottom")

    fig.suptitle("% of ORFs matching TransCode proteomics data")
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(mo):
    mo.md("""
    ## ROC vs TransCode
    """)
    return


@app.cell
def _(GENOME_FASTA, extract_gtf_cds_proteins, ribotie_gtf, trimmed):
    ribotie_proteins = extract_gtf_cds_proteins(ribotie_gtf, GENOME_FASTA)
    ribotish_proteins = extract_gtf_cds_proteins(trimmed, GENOME_FASTA)
    return ribotie_proteins, ribotish_proteins


@app.cell
def _(
    ORFannotate_summary,
    build_roc_dfs,
    phase1_seqs,
    phase2_seqs,
    ribotie_gtf,
    ribotie_proteins,
    ribotish_proteins,
):
    roc_ribotie_p1, roc_ribotish_p1 = build_roc_dfs(
        ribotie_gtf, ORFannotate_summary, ribotie_proteins, ribotish_proteins, phase1_seqs
    )
    roc_ribotie_p2, roc_ribotish_p2 = build_roc_dfs(
        ribotie_gtf, ORFannotate_summary, ribotie_proteins, ribotish_proteins, phase2_seqs
    )
    return roc_ribotie_p1, roc_ribotie_p2, roc_ribotish_p1, roc_ribotish_p2


@app.cell
def _(
    plt,
    roc_auc_score,
    roc_curve,
    roc_ribotie_p1,
    roc_ribotie_p2,
    roc_ribotish_p1,
    roc_ribotish_p2,
):
    _, ax_roc = plt.subplots(figsize=(7, 6))

    for _scores, _labels, _name, _ls in [
        (roc_ribotie_p1["ribotie_score"],  roc_ribotie_p1["label"],  "RiboTIE – Phase 1",   "-"),
        (roc_ribotie_p2["ribotie_score"],  roc_ribotie_p2["label"],  "RiboTIE – Phase 2",   "--"),
        (roc_ribotish_p1["score"],         roc_ribotish_p1["label"], "Ribo-TISH – Phase 1", "-"),
        (roc_ribotish_p2["score"],         roc_ribotish_p2["label"], "Ribo-TISH – Phase 2", "--"),
    ]:
        _fpr, _tpr, _ = roc_curve(_labels, _scores, pos_label=1)
        _auc = roc_auc_score(_labels, _scores)
        ax_roc.plot(_fpr, _tpr, linestyle=_ls, label=f"{_name} (AUC={_auc:.2f})")

    ax_roc.plot([0, 1], [0, 1], "k--")
    ax_roc.set_xlabel("False Positive Rate")
    ax_roc.set_ylabel("True Positive Rate")
    ax_roc.set_title("ROC: recapturing TransCode proteomics ORFs")
    ax_roc.legend()
    plt.tight_layout()
    plt.show()
    return


if __name__ == "__main__":
    app.run()
