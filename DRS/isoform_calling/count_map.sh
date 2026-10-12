#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=count_map
#SBATCH --output=logs/count_map%j.log
#SBATCH --array=1-14

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir count_maps/

config=map_config.txt

sampleid=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
strand=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

python make_count_map.py ${sampleid} ${strand}
