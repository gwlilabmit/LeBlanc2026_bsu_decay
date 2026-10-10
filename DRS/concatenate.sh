#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --job-name=concatenate
#SBATCH --mail-type=END
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/concatenate_%j.log
#SBATCH --array=1-10

mkdir ../alignments
mkdir ../alignments/fastcat_reports/
mkdir ../alignments/concatenated_fastqs/
mkdir ../alignments/split_fastqs

config=concatenate_config.txt

sampleid=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
path_to_fastqs=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
split_prefix=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)

fastcat -f ../alignments/fastcat_reports/${sampleid}_per_file.txt -H -r ../alignments/fastcat_reports/${sampleid}_per_read.txt ${path_to_fastqs} > ../alignments/concatenated_fastqs/${sampleid}_concatenated.fastq

split -l 20000000 -d --additional-suffix=.fastq ../alignments/concatenated_fastqs/${sampleid}_concatenated.fastq ${split_prefix}