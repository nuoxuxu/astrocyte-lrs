#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o sbatch_merged_aligned_flnc_bambu.out
#SBATCH -J sbatch_merged_aligned_flnc_bambu

eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate ./env

bin/run_bambu.R \
/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf \
'nextflow_results/align/minimap2/merged*.minimap.bam' \
/project/rrg-shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa \
192 \
merged_aligned_flnc_bambu.rds