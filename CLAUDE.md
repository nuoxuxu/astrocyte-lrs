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
- Other available subworkflows in this stage: `quality` (GET_QUALITY_METRICS), `ribotie_postanalysis`, `cds_length_distribution`, `vep` (RUN_VEP).

### Key Design Patterns

- **Single-filter strategy**: `min_reads` (default: 5) and `min_n_sample` (default: 2) define the expression filter. This is the only filter level; outputs go directly into their respective `nextflow_results/` subdirectories without a stringency-level subfolder.
- **JSON manifests**: Cross-workflow communication between `main.nf`, `RiboTIE.nf`, and `post_RiboTIE.nf` via JSON files in `nextflow_results/manifests/`.
- **storeDir**: Processes use `storeDir` for persistent output caching (not Nextflow's default `publishDir`).

Process labels control SLURM resource allocation: `short_slurm_job` (1h), `mid_slurm_job` (4h), `long_slurm_job` (24h).

### Scripts

- `bin/` - Executable scripts automatically added to Nextflow PATH. Python scripts use `bin/src/` for shared utilities (GTF parsing, logging).
- `scripts/` - Analysis scripts not called by the pipeline (ad-hoc exploration, SLURM batch templates in `scripts/sbatch/`).

### Environment

Managed via Conda (`environment.yml`) and Apptainer/Singularity containers. Key dependencies: isoseq, samtools, STAR, polars, gffutils, oarfish, and R with IsoformSwitchAnalyzeR.

## Output Structure

Results go to `params.outdir` (default: `nextflow_results`), organized by stage (e.g., `sqanti3/isoseq/sqanti3_filter/`).
