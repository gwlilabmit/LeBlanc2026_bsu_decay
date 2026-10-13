#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=norm_counts
#SBATCH --output=logs/norm_counts%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

echo "python combine_normalize_iso_counts.py"

python combine_normalize_iso_counts.py