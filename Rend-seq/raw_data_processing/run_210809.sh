#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M mtien@mit.edu
#$ -pe whole_nodes 4
#############################################

# removing 3’ adapter
nohup python removeLinker_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7347_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7348_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7349_NA_sequence.fastq
nohup python removeLinker_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7350_NA_sequence.fastq

# aligning to B. subtilis genome
module add bowtie
nohup python align_sequences_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7347_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7348_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7349_NA_sequence-stripped_v2.fastq bsub v1k1 
nohup python align_sequences_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7350_NA_sequence-stripped_v2.fastq bsub v1k1 


# generating wig file
nohup python align_to_wig_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7347_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210809/fastq/210809Li_D21-7347_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7348_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210809/fastq/210809Li_D21-7348_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7349_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210809/fastq/210809Li_D21-7349_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub
nohup python align_to_wig_mzt.py /home/mtien/data/210809/fastq/210809Li_D21-7350_NA_sequence-stripped_v2_bsubtilis_v1k1.align /home/mtien/data/210809/fastq/210809Li_D21-7350_NA_sequence-stripped_v2_bsubtilis_v1k1 bsub