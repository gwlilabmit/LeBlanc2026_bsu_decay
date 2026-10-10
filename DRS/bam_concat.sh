#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --job-name=bam_concat
#SBATCH --mail-type=END
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/bam_concat_%j.log
#SBATCH --array=1-10

mkdir ../alignments/sam_files/
mkdir ../alignments/bam_files/
mkdir ../alignments/flagstat/

config=sam_config.txt
sampleid=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

samtools merge -f ../alignments/bam_files/${sampleid}.bam ../alignments/split_alignments/${sampleid}*.bam

echo 'bam files merged'

samtools sort -o ../alignments/bam_files/${sampleid}.sorted.bam ../alignments/bam_files/${sampleid}.bam

echo 'bam files sorted'

samtools index ../alignments/bam_files/${sampleid}.sorted.bam

echo 'bam files indexed'

samtools view ../alignments/bam_files/${sampleid}.sorted.bam > ../alignments/sam_files/${sampleid}.sorted.sam

echo 'sam files created'

igvtools index ../alignments/sam_files/${sampleid}.sorted.sam

echo 'sam files indexed'

samtools flagstat ../alignments/bam_files/${sampleid}.sorted.bam > ../alignments/flagstat/${sampleid}_flagstat.log

echo 'flagstat done'