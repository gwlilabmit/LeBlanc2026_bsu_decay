#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/full_pipeline_%j.log

#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --job-name=pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leblanch@mit.edu
#SBATCH --output=logs/full_pipeline_%j.log

## quantify all normalization factors in all samples
python count_spike_ins.py
python count_ncRNAs.py
echo 'normalization factors counted'

## normalize ncRNAs by spike-ins
python normalize_ncRNAs.py
echo 'ncRNAs normalized'

## normalization factor coefficients of variation
python normalization_cvs.py
echo 'normalization CVs calculated'

## bulk mRNA quantification
python 190710_bulk_mRNA.py
echo '190710 bulk mRNA quantified'

## cggR fit evaluation
python cggR_fit_eval.py
echo 'cggR fit evaluated'

## differentially stable cggR regions
python cggR_fragment.py
echo 'cggR fragments deconvoluted'

## fit datasets
python wt_fit.py
echo 'wild type fit'

python rnjAdep_fit.py
echo 'rnjA depletion fit'

python rnjAexp_fit.py
echo 'rnjA expression fit'

## filtering
python fits/filtering.py
echo 'datasets filtered for fit quality'

###
echo 'done'