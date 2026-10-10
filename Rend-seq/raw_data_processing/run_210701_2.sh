#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M mtien@mit.edu
#$ -pe whole_nodes 4
#############################################

# removing 3’ adapter
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222008_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222009_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222010_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222011_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222012_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222013_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222014_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222015_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222016_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222017_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222018_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222019_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222020_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222021_NA_sequence.fastq 
nohup python removeLinker_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222022_NA_sequence.fastq 

# aligning to B. subtilis genome
module add bowtie
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222008_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222009_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222010_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222011_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222012_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222013_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222014_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222015_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222016_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222017_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222018_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222019_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222020_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222021_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222022_NA_sequence-stripped_v2.fastq bsub v1k1
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222023_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222024_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222025_NA_sequence-stripped_v2.fastq bsub v1k1


# generating wig file
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222008_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222008_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222009_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222009_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222010_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222010_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222011_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222011_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222012_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222012_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222013_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222013_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222014_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222014_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222015_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222015_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222016_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222016_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222017_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222017_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222018_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222018_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222019_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222019_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222020_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222020_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222021_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222021_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210701/fastq/210701Li_D21-222022_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210701/fastq/210701Li_D21-222022_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
