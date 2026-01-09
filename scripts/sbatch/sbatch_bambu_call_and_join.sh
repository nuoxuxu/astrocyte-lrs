#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o slurm_logs/sbatch_bambu_call_and_join.out
#SBATCH -J sbatch_bambu_call_and_join

eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate ./env

scripts/run_bambu_high_memory.R \
/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf \
'nextflow_results/align/minimap2/merged_IsoSeqX_bc01_5p--IsoSeqX_3p.minimap.bam' \
/project/rrg-shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa \
192 \
bc01.rds