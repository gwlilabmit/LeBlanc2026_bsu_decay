#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=normalize
#SBATCH --output=logs/normalize_%j.log
#SBATCH --array=1-7

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

winsor=$1

DIR=normalized_maps_${winsor}/
mkdir $DIR

config=normalize_config_samples.txt

samplename=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "sbatch normalize_jobs.sh ${samplename} f ${winsor}"

sbatch normalize_jobs.sh ${samplename} f ${winsor}

echo "sbatch normalize_jobs.sh ${samplename} r ${winsor}"

sbatch normalize_jobs.sh ${samplename} r ${winsor}
