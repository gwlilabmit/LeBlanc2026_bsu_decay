#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M mtien@mit.edu
#$ -pe whole_nodes 4
#############################################

# removing 3’ adapter
FILES=/home/mtien/data/190710/fastq/190710Li_D19-7*.fastq
##FILES=/home/mtien/data/190710/fastq/190710Li_D19-7618_NA_sequence.fastq
for f in $FILES
do
        nohup python removeLinker_mzt.py $f 
done

# aligning to B. subtilis genome
FILES=/home/mtien/data/190710/fastq/190710Li_D19-7*stripped_v2.fastq
module add bowtie

##FILES=/home/mtien/data/190710/fastq/190710Li_D19-7618_NA_sequence-stripped_v2.fastq
for f in $FILES
do
        nohup python align_sequences_mzt.py $f bsub v1k1
done

# generating wig file
FILES=/home/mtien/data/190710/fastq/190710Li_D19-7*.align
##FILES=/home/mtien/data/190710/fastq/190710Li_D19-7618_NA_sequence-stripped_v2_bsubtilis_v1k1.align
for f in $FILES
do
        nohup python align_to_wig_mzt.py $f $f bsub
done
