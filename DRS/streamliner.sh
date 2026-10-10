#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH --job-name=streamline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/streamline_%j.log
#SBATCH --array=1-8

source ~/miniconda3/etc/profile.d/conda.sh
conda activate isoforms

mkdir streamlined_sam_files/

config=sam_config.txt
sampleid=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

python sam_streamliner.py --infile ~/data/nanopore/250603Li/alignments/sam_files/${sampleid}.sorted.sam --outfile streamlined_sam_files/${sampleid}.txt
