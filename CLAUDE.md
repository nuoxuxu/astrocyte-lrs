# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**astrocyte-lrs** is a Nextflow DSL2 pipeline for analyzing long-read RNA sequencing (Iso-Seq) data from astrocyte cells. It performs transcriptome analysis including read preprocessing, transcript discovery, expression quantification, quality assessment, and protein-coding region prediction.

## Running the Pipeline

```bash
# Main pipeline (HPC with SLURM)
nextflow run main.nf -profile trillium

# RiboTIE training (requires GPU)
nextflow run RiboTIE.nf -profile trillium_gpu

# Local execution
nextflow run main.nf -profile local
```

Profiles: `trillium` (local HPC), `trillium_gpu` (GPU jobs), `narval` (Narval cluster), `local` (single machine). All parameters are in `nextflow.config`.

## Architecture

### Pipeline Stages (main.nf)

The pipeline is a linear chain of subworkflows, each in `subworkflows/<name>/main.nf`:

```
PREPROCESSING → ISOSEQ → RUN_OARFISH → SQANTI_AND_FILTER_BY_EXP → RUN_ORFANAGE → PREPARE_RIBOTIE
                                                                         ↓
                                                          GET_QUALITY_METRICS (PhyloCSF++, Pfam)
                                                          ISOFORMSWITCH (R-based)
                                                          RIBOTIE_VISUALIZATION
```

1. **preprocessing** - PacBio demux (skera → lima → refine) producing FLNC BAMs
2. **isoseq** - Transcript clustering and genome alignment (cluster2 → pbmm2 → collapse)
3. **oarfish** - Transcript quantification via minimap2 + oarfish
4. **sqanti** - SQANTI3 QC/filtering + expression-based filtering with dual configs
5. **orfanage** - ORF annotation and protein extraction
6. **riboseq** - Ribo-seq alignment and RiboTIE database preparation
7. **quality** - PhyloCSF++ conservation scores and Pfam domain scanning
8. **IsoformSwitchAnalyzeR** - Differential isoform usage analysis (R)
9. **visualization** - RiboTIE publication figures

### Key Design Patterns

- **Single-filter strategy**: `min_reads` (default: 5) and `min_n_sample` (default: 2) define the expression filter applied as `mid_stringency`. Outputs are tagged with `param_set_name = "mid_stringency"` as a tuple element.
- **JSON manifests**: Cross-workflow communication between `main.nf`, `post_RiboTIE,nf` and `RiboTIE.nf` via JSON files in `nextflow_results/manifests/`.
- **storeDir**: Processes use `storeDir` for persistent output caching (not Nextflow's default `publishDir`).
- **Channel tuples**: Data flows as `[param_set_name, file]` tuples for tracking filter config provenance.

Process labels control SLURM resource allocation: `short_slurm_job` (1h), `mid_slurm_job` (4h), `long_slurm_job` (24h).

### Scripts

- `bin/` - Executable scripts automatically added to Nextflow PATH. Python scripts use `bin/src/` for shared utilities (GTF parsing, logging).
- `scripts/` - Analysis scripts not called by the pipeline (ad-hoc exploration, SLURM batch templates in `scripts/sbatch/`).

### Environment

Managed via Conda (`environment.yml`) and Apptainer/Singularity containers. Key dependencies: isoseq, samtools, STAR, polars, gffutils, oarfish, and R with IsoformSwitchAnalyzeR.

## Output Structure

Results go to `params.outdir` (default: `nextflow_results`), organized by stage (e.g., `sqanti3/isoseq/sqanti3_filter/mid_stringency/`).
