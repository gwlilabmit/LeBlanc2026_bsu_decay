#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --nice=100000
#SBATCH --exclude=c[5-22]
#SBATCH --job-name=subnorm
#SBATCH --output=logs/subnorm_%j.log
#SBATCH --array=1-16

samplename=$1
strand=$2
winsor=$3

if [ $strand == 'f' ]
then
    config=normalize_config_sections_f.txt
elif [ $strand == 'r' ]
then
    config=normalize_config_sections_r.txt
else
    echo 'must enter valid strand'
fi

sectionname=$(awk -F, -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

echo "cross_normalize.py ${samplename}_${strand} ${sectionname} ${winsor}"

python cross_normalize.py ${samplename}_${strand} ${sectionname} ${winsor}

wait
echo "done"