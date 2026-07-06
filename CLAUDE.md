# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**astrocyte-lrs** is a Nextflow DSL2 pipeline for analyzing long-read RNA sequencing (Iso-Seq) data from astrocyte cells. It performs transcriptome analysis including read preprocessing, transcript discovery, expression quantification, quality assessment, and protein-coding region prediction.

## Running the Pipeline

The pipeline runs in three sequential entrypoints, connected by JSON manifests:

```bash
# 1. Main pipeline: preprocessing → RiboTIE database prep (HPC with SLURM)
nextflow run main.nf -profile trillium

# 2. RiboTIE training/inference (requires GPU; consumes manifests from step 1)
nextflow run RiboTIE.nf -profile trillium_gpu

# 3. Post-RiboTIE downstream analysis (IsoformSwitch, aim_2, quality, summary table)
nextflow run post_RiboTIE.nf -profile trillium

# Local execution
nextflow run main.nf -profile local
```

Profiles: `trillium` (local HPC), `trillium_gpu` (GPU jobs), `narval` (Narval cluster), `local` (single machine). All parameters are in `nextflow.config`.

## Architecture

Subworkflows live in `subworkflows/local/<name>/main.nf`.

### Stage 1 — `main.nf` (transcript discovery → RiboTIE prep)

```
PREPROCESSING → ISOSEQ → RUN_OARFISH → SQANTI → FILTER_BY_EXPRESSION → RUN_ORFANAGE → PREPARE_RIBOTIE
```

1. **preprocessing** - PacBio demux (skera → lima → refine) producing FLNC BAMs
2. **isoseq** - Transcript clustering and genome alignment (cluster2 → pbmm2 → collapse)
3. **oarfish** - Transcript quantification via minimap2 + oarfish
4. **sqanti** - SQANTI3 QC/filtering
5. **filter_by_expression** - Expression-based filtering of transcripts
6. **orfanage** - ORF annotation and protein extraction
7. **prepare_ribotie** - Ribo-seq alignment and RiboTIE database (h5) preparation; writes JSON manifests to `nextflow_results/manifests/`

### Stage 2 — `RiboTIE.nf` (GPU)

RiboTIE training/inference, driven by the manifests written by stage 1.

### Stage 3 — `post_RiboTIE.nf` (downstream analysis)

**Inputs come from the `from_collaborator/` folder, not from a stage-2 output in this repo.** The RiboTIE ORF predictions were produced by a collaborator (the code that derived them is not in this repo). Key files:

- `from_collaborator/ribotie_cpm1_3sample.csv` — per-ORF RiboTIE predictions (ORF_id, transcript_id, scores, start/stop codons, coordinates).
- `from_collaborator/filtered_output.gtf` — GTF of the predicted ORFs (CDS/exon/start_codon/stop_codon features, keyed by `ORF_id`, where `transcript_id == ORF_id`).
- `from_collaborator/targetable_asos_passQC.csv` — ASO-targetable ORFs filtered from `summary_table.tsv` for feasibility of antisense oligonucleotide perturbation. Columns include: `gene_symbol`, `orf_id`, `transcript`, `ORF_type`, `gene_type`, `region_type`, `chrom`/`start`/`end`/`strand`, `seq` (20-nt target region), `aso_antisense` (the ASO sequence), GC/homopolymer/CpG/self-complement QC metrics, `n_sibling_hits`/`n_off` (off-target counts), and `qc_pass`/`qc_flags`. This is the final prioritized list for experimental ASO perturbation to test for phenotypic changes in cultured astrocytes.

The RiboTIE files were derived from the stage-1 `prepare_ribotie` output run through RiboTIE by the collaborator, so treat them as the authoritative ORF set for downstream analysis. Stage 3 also reuses stage-1 outputs (e.g. `final_*` under `nextflow_results/sqanti3/...`).

Downstream flow:

```
supplement_collaborator_gtf → ISOFORMSWITCH (R) → AIM_2 (sQTL/coloc)
                            → FILTER_RIBOTIE → SUMMARY_TABLE
```

- **IsoformSwitchAnalyzeR** - Differential isoform usage analysis (R); the `ALL2` run feeds the summary table. In `isoformFeatures.csv`, `gene_switch_q_value` is the minimum `isoform_switch_q_value` across all isoforms of a gene — it is not an independent gene-level test.
- **aim_2** - Novel coding junctions and leafcutter/novel-junction sQTL–coloc matching
- **filter_ribotie** - Filters RiboTIE ORFs/proteins for downstream scoring
- **summary_table** - Per-ORF evidence table (`nextflow_results/summary_table/summary_table.tsv`); combines CPAT, Pfam, PhyloCSF++, GENCODE/Study2 novelty, IsoformSwitch metrics, sQTL coloc, and Ribo-seq concordance. Column meanings are documented in `docs/summary_table_columns.md`. Built by `bin/make_summary_table.py`.
- Other available subworkflows in this stage: `ribotie_postanalysis`, `cds_length_distribution`, `vep` (RUN_VEP).

### `quality.nf` — ORF quality metrics across parameter sets

A standalone entrypoint (`nextflow run quality.nf -profile trillium`) that adds quality
annotations to RiboTIE-predicted ORFs for each Iso-Seq parameter set (minlen, no_minlen).
It pairs `ribotie_training_outputs_{name}.json` (RiboTIE predictions) with the matching
`main_pipeline_outputs_{name}.json` (stage-1 reference data) by name. The `gencode`
parameter set has no `main_pipeline_outputs_gencode.json` and is skipped in workflows that
require Iso-Seq auxiliary data.

Current workflows exported from `subworkflows/local/quality/main.nf`:
- **`GET_QUALITY_METRICS`** - PhyloCSF++ conservation annotation of novel ORF GTFs
- **`LABEL_ORF_TYPE_GENCODE`** - Adds `ORF_type_ORFanage` and `ORF_type_RiboTIE` to the RiboTIE
  merged CSV. Output: `nextflow_results/quality/{name}_orf_type_gencode.tsv`.
  Also produces `{name}_final_quality_metrics.tsv` with both `ORF_type_ORFanage` and
  `ORF_type_RiboTIE` inserted immediately before `ORBLv` in the riboseq quality metrics table.

#### ORF type columns — semantics (`bin/label_orf_type_gencode.py`)

Both columns classify an ORF against the same reference (GENCODE canonical CDS projected onto the Iso-Seq transcript's exon structure), but differ in which ORF is the **query**:

| Column | Query | Reference | Purpose |
|--------|-------|-----------|---------|
| `ORF_type_ORFanage` | ORFanage-predicted CDS for this transcript | GENCODE canonical CDS (Tier-1 only: `orfanage_template`) | Reflects the ORF type label seen by RiboTIE during training |
| `ORF_type_RiboTIE` | RiboTIE-predicted ORF (from `TIS_pos`/`TTS_pos`) | GENCODE canonical CDS (three-tier lookup) | How the RiboTIE prediction compares to GENCODE |

`ORF_type_ORFanage` uses only Tier-1 because that is exactly what RiboTIE was trained on — ORFanage derives its CDS from the `orfanage_template` ENST, and that template is what was used to generate training labels. `ORF_type_RiboTIE` uses all three tiers to give the best-available GENCODE reference for every transcript, including those ORFanage did not template.

#### In-frame vs ncORF categories

ORF types are divided into two groups:

**In-frame** — the ORF shares the reading frame of a GENCODE canonical CDS (same TIS or same TTS, so codon phase is preserved):

| Category | Relationship to GENCODE canonical |
|----------|----------------------------------|
| `annotated CDS` | Identical start and stop |
| `N-terminal extension` | Earlier start, same stop |
| `C-terminal extension` | Same start, later stop |
| `N-terminal truncation` | Later start, same stop |
| `C-terminal truncation` | Same start, earlier stop |

**ncORF** — everything else:
- `uORF`, `uoORF`, `dORF`, `doORF`, `intORF` — different frame or entirely outside the canonical CDS
- `lncRNA-ORF` — the resolved ENST is a lncRNA; no canonical CDS exists, so there is no frame to be in. Classified as ncORF even though the ORF occupies a transcribed region.
- `varRNA-ORF` — same rationale; resolved ENST is another non-coding biotype
- `no_template_*` — no GENCODE CDS template could be resolved at all

`plot_orf_type_bars.py` (called by the `plot_orf_type_bars` process) shows only in-frame ORFs; percentages are computed over the in-frame total per column, not over all ORFs.

#### Non-coding fallback: lncRNA-ORF vs varRNA-ORF

When no GENCODE canonical CDS can be projected (missing template, non-coding ENST, or TIS outside all exons), the ORF is labelled based on the resolved ENST's `transcript_type` in GENCODE:
- `lncRNA-ORF` if `transcript_type == 'lncRNA'`
- `varRNA-ORF` for all other non-coding biotypes (or when no ENST is resolved at all)

This mirrors RiboTIE's own classification logic from training.

#### Three-tier canonical CDS lookup (for `ORF_type_RiboTIE`)

For each Iso-Seq transcript, the GENCODE canonical CDS is resolved via a prioritised lookup before running `project_gencode_cds` / `classify_orf_type`:

| Tier | Source | Condition |
|------|--------|-----------|
| 1 | `orfanage_template` attribute in `orfanage.gtf` | Transcript appears in the plain (non-numbered-exons) orfanage GTF with an `orfanage_template` attribute |
| 2 | `associated_transcript` from SQANTI3 `final_classification.parquet` | `associated_transcript` is not `'novel'` |
| 3 | MANE_Select transcript for `associated_gene` in GENCODE | `associated_transcript == 'novel'` and the gene has a MANE_Select transcript |

ENSEMBL version suffixes are stripped before all lookups (e.g. `ENST00000361390.2` → `ENST00000361390`) so that version mismatches between sources do not cause false misses.

If the resolved ENST has no CDS annotation in GENCODE, the ORF receives a `no_template_*` label instead of being classified. All special-case labels:

| Label | Condition |
|-------|-----------|
| `no_template_fusion` | `associated_gene` matches `ENSG…_ENSG…` pattern (fusion transcript) |
| `no_template_novel_gene` | `associated_gene` starts with `novelGene_` |
| `no_template_unknown_gene` | transcript not found in `final_classification` |
| `no_template_no_MANE_Select` | `associated_transcript == 'novel'` but no MANE_Select exists for the gene |
| `no_template_orfanage_template_non_coding` | Tier-1 ENST found but has no CDS in GENCODE |
| `no_template_associated_transcript_non_coding` | Tier-2 ENST found but has no CDS in GENCODE |
| `no_template_associated_MANE_Select_non_coding` | Tier-3 MANE_Select ENST found but has no CDS in GENCODE |

Plain orfanage GTF paths are hardcoded in `nextflow.config` as `params.orfanage_gtf_no_minlen` and `params.orfanage_gtf_minlen`.

### Key Design Patterns

- **Single-filter strategy**: `min_reads` (default: 5) and `min_n_sample` (default: 2) define the expression filter. This is the only filter level; outputs go directly into their respective `nextflow_results/` subdirectories without a stringency-level subfolder.
- **JSON manifests**: Cross-workflow communication between `main.nf`, `RiboTIE.nf`, `post_RiboTIE.nf`, and `quality.nf` via JSON files in `nextflow_results/manifests/`. Two manifest types exist per parameter set (minlen, no_minlen): `main_pipeline_outputs_{name}.json` (stage-1 outputs: orfanage_gtf, final_classification, final_fasta, final_expression) and `ribotie_training_outputs_{name}.json` (RiboTIE predictions: ribotie_merged_csv/gtf and novel/redundant variants). A `ribotie_training_outputs_gencode.json` also exists for the GENCODE-reference RiboTIE run but has no matching `main_pipeline_outputs_gencode.json`.
- **storeDir**: Processes use `storeDir` for persistent output caching (not Nextflow's default `publishDir`).

Process labels control SLURM resource allocation: `short_slurm_job` (1h), `mid_slurm_job` (4h), `long_slurm_job` (24h).

### Scripts

- `bin/` - Executable scripts automatically added to Nextflow PATH. Python scripts use `bin/src/` for shared utilities (GTF parsing, logging).
- `scripts/` - Analysis scripts not called by the pipeline (ad-hoc exploration, SLURM batch templates in `scripts/sbatch/`).

### Environment

Managed via Conda (`environment.yml`) and Apptainer/Singularity containers. Key dependencies: isoseq, samtools, STAR, polars, gffutils, oarfish, and R with IsoformSwitchAnalyzeR.

## Output Structure

Results go to `params.outdir` (default: `nextflow_results`), organized by stage (e.g., `sqanti3/isoseq/sqanti3_filter/`).
