# astrocyte-lrs

Nextflow DSL2 pipeline for analyzing long-read RNA sequencing (Iso-Seq) and Ribo-seq data from astrocyte cells. Performs transcriptome analysis including read preprocessing, transcript discovery, expression quantification, quality assessment, and protein-coding region prediction.

## Data

### Long-read RNA-seq
- PacBio Kinnex bulk RNA-seq on Revio platform (4-plex)
- Input: HiFi reads BAM files (`data/long_read/pacbio/`)

### Short-read RNA-seq
- Library prepared using Takara stranded mRNA kit, sequenced on NovaSeq X platform
- Input: paired-end FASTQ files (`data/short_read/`)

### Ribo-seq
- Reads unmapped to contaminants used as input
- Input: `data/ribo_seq/`

### Reference data

Download external reference datasets:

```bash
# TransCODE phase 1 data
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01369-0/MediaObjects/41587_2022_1369_MOESM2_ESM.xlsx -O data/41587_2022_1369_MOESM2_ESM.xlsx

# TransCODE phase 2 data
wget https://www.biorxiv.org/content/biorxiv/early/2025/07/07/2025.07.03.662928/DC1/embed/media-1.xlsx -O data/media-1.xlsx
```

## Installation

### Conda environment

```bash
conda env create -f environment.yml
conda activate ./env
```

### SQANTI3

```bash
wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.5.4/SQANTI3_v5.5.4.zip
mkdir sqanti3
unzip SQANTI3_v5.5.4.zip -d sqanti3
```

Key dependencies: `isoseq`, `samtools`, `STAR`, `polars`, `gffutils`, `oarfish`, and R with IsoformSwitchAnalyzeR. Processes run inside Apptainer/Singularity containers; container images are cached at `NXF_SINGULARITY_CACHEDIR`.

## Configuration

All parameters are defined in `nextflow.config`. Key parameters:

| Parameter | Description |
|---|---|
| `hifi_reads_bam` | Glob pattern for PacBio HiFi BAM files |
| `short_read_fastqs` | Glob pattern for short-read paired FASTQ files |
| `riboseq_unmapped_to_contaminants` | Glob pattern for Ribo-seq unmapped reads |
| `ref_genome_fasta` | GRCh38 reference genome FASTA |
| `annotation_gtf` | GENCODE annotation GTF (v47) |
| `biosamples_csv` | Sample metadata for demultiplexing |
| `filter_configs` | List of `[name, min_reads, min_samples]` filter stringency sets |
| `outdir` | Output directory (default: `nextflow_results`) |

### Filter stringency sets

Three parallel filter configurations are run through the pipeline:

| Name | Min reads | Min samples |
|---|---|---|
| `low_stringency` | 1 | 2 |
| `mid_stringency` | 5 | 2 |
| `high_stringency` | 5 | 3 |

Downstream steps use `mid_stringency` results by default (e.g., IsoformSwitchAnalyzeR, post-RiboTIE analysis).

## Running the Pipeline

The analysis is split across three Nextflow workflows that run sequentially:

```bash
# Step 1: Main pipeline — preprocessing through RiboTIE database preparation
nextflow run main.nf -profile trillium

# Step 2: RiboTIE training — requires GPU
nextflow run RiboTIE.nf -profile trillium_gpu

# Step 3: Post-RiboTIE — quality metrics, isoform switching, visualization
nextflow run post_RiboTIE.nf -profile trillium
```

Workflows communicate via JSON manifests written to `nextflow_results/manifests/`:
- `ribotie_training_inputs.json` — inputs for `RiboTIE.nf`
- `ribotie_training_outputs.json` — RiboTIE result paths for `post_RiboTIE.nf`
- `main_pipeline_outputs.json` — filtered transcriptome paths for `post_RiboTIE.nf`

### Execution profiles

| Profile | Description |
|---|---|
| `trillium` | Local HPC with SLURM (Trillium cluster) |
| `trillium_gpu` | GPU jobs on Trillium (for RiboTIE training) |
| `narval` | Narval cluster (Alliance Canada) |
| `local` | Single machine, no SLURM |

### Local execution

```bash
nextflow run main.nf -profile local
```

## Pipeline Stages

### main.nf

```
PREPROCESSING → ISOSEQ → RUN_OARFISH → SQANTI → FILTER_BY_EXPRESSION → RUN_ORFANAGE → PREPARE_RIBOTIE
```

1. **PREPROCESSING** — PacBio demultiplexing: `skera` (Kinnex adapter removal) → `lima` (primer removal) → `isoseq refine` (FLNC BAM generation)
2. **ISOSEQ** — Transcript clustering and genome alignment: `cluster2` → `pbmm2` → `collapse`
3. **RUN_OARFISH** — Transcript quantification: `minimap2` alignment + `oarfish` EM quantification
4. **SQANTI** — SQANTI3 QC and filtering using short-read junction support + expression-based filtering
5. **FILTER_BY_EXPRESSION** — Applies the three stringency filter configs in parallel; outputs tagged `[param_set_name, file]`
6. **RUN_ORFANAGE** — ORF annotation and protein sequence extraction with ORFanage
7. **PREPARE_RIBOTIE** — Ribo-seq alignment (STAR) and RiboTIE database preparation

### RiboTIE.nf

Trains and runs the RiboTIE deep-learning model on prepared databases. Reads inputs from `ribotie_training_inputs.json`. Requires a GPU node (`trillium_gpu` profile). Outputs redundant, novel, and filtered ORF predictions per stringency set.

### post_RiboTIE.nf

```
GET_QUALITY_METRICS    — PhyloCSF++ conservation scores + Pfam domain scanning
ISOFORMSWITCH          — Differential isoform usage (IsoformSwitchAnalyzeR, R)
RIBOTIE_POSTANALYSIS   — RiboTIE result post-processing and visualization
```

## Output Structure

Results are written to `params.outdir` (default: `/scratch/nxu/astrocytes/nextflow_results`), organized by stage and stringency level:

```
nextflow_results/
├── manifests/                         # JSON cross-workflow manifests
├── preprocessing/                     # FLNC BAMs per sample
├── isoseq/                            # Collapsed transcript GTF + BAM
├── oarfish/                           # Quantification results
├── sqanti3/isoseq/sqanti3_filter/
│   ├── low_stringency/
│   ├── mid_stringency/
│   └── high_stringency/
├── orf_prediction/orfanage/           # ORF GTF + protein FASTA per stringency
├── prepare_ribotie/                   # STAR BAMs + RiboTIE .h5 databases
└── ribotie/                           # RiboTIE predictions per stringency
```

Processes use `storeDir` for persistent caching; re-running skips completed processes automatically.

## Scripts

- `bin/` — Executable scripts on the Nextflow PATH. Shared utilities (GTF parsing, logging) are in `bin/src/`.
- `scripts/` — Ad-hoc analysis scripts and SLURM batch templates (`scripts/sbatch/`). Not called by the pipeline.
  - `scripts/deduplicate_cds_transcripts.py` — Deduplicates RiboTIE transcripts that share identical CDS coordinates; writes `transcript_id` → `transcript_id_base` mapping CSVs to `export/local/`.
