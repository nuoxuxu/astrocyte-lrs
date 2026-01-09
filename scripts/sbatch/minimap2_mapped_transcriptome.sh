#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o slurm_logs/minimap2_mapped_transcriptome.out
#SBATCH -J minimap2_mapped_transcriptome

eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate ./env

minimap2 --eqx -N 100 -ax map-hifi \
    -t 192 proc/joint_transcriptome.fasta \
    nextflow_results/prepare/convert_ubam_to_fastqz/merged_IsoSeqX_bc01_5p--IsoSeqX_3p.fastq.gz