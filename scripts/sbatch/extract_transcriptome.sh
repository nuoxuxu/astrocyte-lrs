#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o slurm_logs/extract_transcriptome.out
#SBATCH -J extract_transcriptome

eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate ./env

gffread -w proc/joint_transcriptome.fasta -g /project/rrg-shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa nextflow_results/sqanti3/isoseq/sqanti3_filter/default.filtered.gtf