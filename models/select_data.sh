#!/bin/bash


#! Give your job a name
#SBATCH -J tensor 
#! How much memory do you need?
#SBATCH --mem=300G
#! How much wallclock time will be required?
#SBATCH --time=1-0
##SBATCH --no-requeue
#SBATCH -p epyc
source ~/.bashrc
conda activate R4_2

Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R \
--ichorcna_tf_cutoff_lower 0  \
--ichorcna_tf_cutoff_upper 1  \
--outdir .  