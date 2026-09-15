#!/bin/sh
#SBATCH --job-name=fragpipe_test
#SBATCH --output=fragpipe_test.out
#SBATCH --error=fragpipe_test.err
#SBATCH --time=4:00:00
#SBATCH -n 1
#SBATCH -N 1

export LC_ALL=C.UTF-8
export LANG=C.UTF-8

module load apptainer

apptainer run \
	--home /scratch/nxu/fragpipe_home:/home/nxu \
	-B /home/nxu/tools:/home/nxu/tools \
	-B /home/nxu/SCRATCH:/home/nxu/SCRATCH \
	-B /scratch/nxu \
	-B /cvmfs/soft.computecanada.ca \
	-B /scratch/nxu/fragpipe_writable/cache:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/cache \
	-B /scratch/nxu/fragpipe_writable/jobs:/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/jobs \
	$NXF_SINGULARITY_CACHEDIR/fragpipe_latest.sif \
	/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin/fragpipe \
	--headless \
	--workflow /scratch/nxu/astrocytes/assets/diaPASEF.workflow \
	--manifest /scratch/nxu/astrocytes/assets/fragpipe_manifest.fp-manifest \
	--workdir /scratch/nxu/astrocytes/results/proteomics \
	--config-tools-folder /home/nxu/tools/fragpipe_tools \
	--config-python /bin/python3