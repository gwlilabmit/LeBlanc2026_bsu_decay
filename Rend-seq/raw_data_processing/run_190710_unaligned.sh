#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ -M mtien@mit.edu
#$ -pe whole_nodes 4
#############################################

# aligning to B. subtilis genome
FILES=/home/mtien/data/190710/fastq/unalign_v1k1/190710Li_D19-7*
module add bowtie
for f in $FILES
do
        nohup python align_sequences_mzt.py $f ivt_spike v2k1
done

# generating wig file
FILES=/home/mtien/data/190710/fastq/unalign_v1k1/190710Li_D19-7*.align
for f in $FILES
do
        nohup python alignIVT_to_wig_mzt.py $f $f ivt_spike
done
