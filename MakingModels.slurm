#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marlo072@umn.edu
#SBATCH --job-name SimModels
#SBATCH -o ./Jobs/%j.SimModels.out
#SBATCH -e ./Jobs/%j.SimModels.err
#SBATCH -p msismall,msilarge

set -e

cd $SLURM_SUBMIT_DIR
conda activate snakeprac
snakemake -s HowMakeModels.smk -c 6\
	--profile slurm.go_wags\
	--use-conda \
	--rerun-incomplete \
