#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o slurm_logs/oarfish_test.out
#SBATCH -J oarfish_test

eval "$(conda shell.bash hook)"
source $CONDA_PREFIX/etc/profile.d/mamba.sh
mamba activate ./env

oarfish --threads 192 \
        --filter-group no-filters \
        --model-coverage \
        --alignments {input} \
        --output {params.output_path}
