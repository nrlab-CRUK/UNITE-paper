

sbatch --time=0-1 --mem=8G --partition="epyc" \
	--wrap="source ~/.bashrc;conda activate R4_3; Rscript /home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_and_cnn.R --which_tf_strat '[0, 0.03]'"

sbatch --time=0-1 --mem=8G --partition="epyc" \
	--wrap="source ~/.bashrc;conda activate R4_3; Rscript /home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_and_cnn.R --which_tf_strat '(0.03, 0.1]'"

sbatch --time=0-1 --mem=8G --partition="epyc" \
	--wrap="source ~/.bashrc;conda activate R4_3; Rscript /home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_and_cnn.R --which_tf_strat '(0.1, 1]'"

sbatch --time=0-1 --mem=8G --partition="epyc" \
	--wrap="source ~/.bashrc;conda activate R4_3; Rscript /home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_and_cnn.R --which_tf_strat 'all'"



