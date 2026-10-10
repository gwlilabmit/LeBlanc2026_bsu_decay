#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=align
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/align_%j.log
#SBATCH --array=1-47

mkdir ../alignments/
mkdir ../alignments/split_alignments/

config=align_config.txt
fastqid=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

minimap2 -a -x map-ont -k14 --MD --secondary=no ../../../genomes/NC_000964.3/B_subtilis_168_NC_000964.3.fa ../alignments/split_fastqs/${fastqid}.fastq > ../alignments/split_alignments/${fastqid}.sam

echo ${fastqid}

samtools view -b -h ../alignments/split_alignments/${fastqid}.sam > ../alignments/split_alignments/${fastqid}.bam