#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --job-name=quant
#SBATCH --output=logs/quant%j.log
#SBATCH --array=1-7

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir isoform_counts_by_timepoint

config=sample_config.txt

timepoint=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "python isoform_quant.py ${timepoint}"

python isoform_quant.py ${timepoint}