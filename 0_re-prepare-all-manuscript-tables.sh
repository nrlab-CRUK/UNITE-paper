
# !/bin/bash
sbatch --time=0-3 --mem=32G --partition="epyc" -J "latex-tables" --wrap="source ~/.bashrc;cd /home/nrlab/wang04/ulyses/; git checkout cuDF-removeChr19-xgb;conda activate R4_3;Rscript /home/nrlab/wang04/ulyses/0_re-prepare-all-manuscript-tables.r"

