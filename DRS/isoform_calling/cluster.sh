#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=cluster
#SBATCH --output=logs/cluster%j.log
#SBATCH --array=1-7

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir isoform_clusters
mkdir dendrograms

config=sample_config.txt

timepoint=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
z_thresh=30

echo "python cluster_isoforms.py ${timepoint} ${z_thresh}"

python cluster_isoforms.py ${timepoint} ${z_thresh} 100

python cluster_isoforms.py ${timepoint} ${z_thresh} 90