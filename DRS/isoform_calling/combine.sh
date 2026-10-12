#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=combine
#SBATCH --output=logs/combine%j.log
#SBATCH --array=1-7

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir combined_isoforms/

config=sample_config.txt
timepoint=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "python combine_winsor.py ${timepoint}"

python combine_winsor.py ${timepoint}
