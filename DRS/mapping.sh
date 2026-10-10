#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH --job-name=mapping
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/mapping_%j.log
#SBATCH --array=1-8

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir mRNA_overlap/

config=sam_config.txt
sampleid=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "mRNA_mapping_reads.py --infile streamlined_sam_files/${sampleid}.txt --outfile mRNA_overlap/${sampleid}.txt"

python mRNA_mapping_reads.py --infile streamlined_sam_files/${sampleid}.txt --outfile mRNA_overlap/${sampleid}.txt
