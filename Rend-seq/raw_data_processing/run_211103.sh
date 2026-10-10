#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M mtien@mit.edu
#$ -pe whole_nodes 4
#############################################

# removing 3’ adapter
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9974_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10011_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10012_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10013_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9978_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10014_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10015_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9981_NA_sequence.fastq 

# aligning to B. subtilis genome
module add bowtie
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9974_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10011_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10012_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10013_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9978_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10014_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10015_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9981_NA_sequence-stripped_v2.fastq bsub v1k1

# generating wig file
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9974_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-9974_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10011_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-10011_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10012_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-10012_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10013_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-10013_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9978_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-9978_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10014_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-10014_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-10015_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-10015_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/211103/fastq/211103Li_D21-9981_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/211103/fastq/211103Li_D21-9981_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub