#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=filter
#SBATCH --output=logs/filter%j.log
#SBATCH --array=1-7

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir filtered_isoforms/

config=sample_config.txt

timepoint=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
z_thresh=30

echo "python filter_isoforms.py ${timepoint} ${z_thresh}"

python filter_isoforms.py ${timepoint} ${z_thresh}