#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --job-name=pipeline
#SBATCH --mail-type=END
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/pipeline_%j.log

source ~/miniconda3/etc/profile.d/conda.sh
conda activate nanopore

# manual step: edit concat_config.txt to match your files
# manual step: set number of array jobs to number of samples in concatenate.sh, bam_concat.sh, streamliner.sh, and mapping.sh

sbatch --wait concatenate.sh

cd ../alignments/split_fastqs
dir > temp.txt
cd ../../scripts
python align_config_maker.py ../alignments/split_fastqs/temp.txt align_config.txt
rm ../alignments/split_fastqs/temp.txt

# # manual step: set number of array jobs in align.sh

sbatch --wait align.sh

cd ../
dir > temp.txt
cd scripts
python sam_config_maker.py ../temp.txt sam_config.txt
rm ../temp.txt

sbatch --wait bam_concat.sh

sbatch --wait streamliner.sh

sbatch --wait mapping.sh

echo 'complete'
