#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --job-name=genes
#SBATCH --output=logs/genes_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir gene_matching/

echo "python gene_matching_all_timepoints.py"

python gene_matching_all_timepoints.py