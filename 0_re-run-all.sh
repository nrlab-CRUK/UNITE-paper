
# !/bin/bash
cd /home/nrlab/wang04/ulyses
git checkout logistic_reg

# plot data split bar plot in Fig. 5a

  sbatch --time=0-2 --mem=16G --partition="epyc" -J "paper-run-datasplit" --wrap="source ~/.bashrc;conda activate R4_3;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg;Rscript /home/nrlab/wang04/ulyses/0_data_split/plot_data_split.r "

#prepare the stat source data for journal
sbatch --time=0-2 --mem=16G --partition="epyc" -J "paper-run-stat-source" --dependency=$(squeue --noheader --format %i --name paper-run-datasplit) --wrap="source ~/.bashrc;conda activate R4_3;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg; Rscript /home/nrlab/wang04/ulyses/6_statistics_source_data/prepare_statistics_raw_data.r "

# direct compare xgb and cnn
  sbatch --time=0-1 --mem=8G --partition="epyc" -J "paper-run-comp-xgb-cnn" --dependency=$(squeue --noheader --format %i --name paper-run-stat-source) --wrap="source ~/.bashrc;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg;bash /home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_cnn.sh "

  sbatch --time=0-1 --mem=8G --partition="epyc" -J "paper-run-comp-xgb-cnn2" --dependency=$(squeue --noheader --format %i --name paper-run-comp-xgb-cnn) --wrap="source ~/.bashrc;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg;conda activate R4_3;Rscript /home/nrlab/wang04/ulyses/1_nested_cv_plots/combine_main_performance_figs.R "

# venn
  sbatch --time=0-1 --mem=8G --partition="epyc" -J "paper-run-venn" --dependency=$(squeue --noheader --format %i --name paper-run-comp-xgb-cnn2) --wrap="source ~/.bashrc;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg;conda activate R4_3;Rscript /home/nrlab/wang04/ulyses/3_Venn/plot_venn_diagram.R"

# indtest scores
  sbatch --time=0-1 --mem=8G --partition="epyc" -J "paper-run-indtest" --dependency=$(squeue --noheader --format %i --name paper-run-venn) --wrap="source ~/.bashrc;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg;bash /home/nrlab/wang04/ulyses/4_ROC/plot_indtest_roc.sh"


#wait 2 mins
sleep 120
# add commit and push to github

  sbatch --time=0-1 --mem=1G --partition="epyc" -J "paper-commit" --dependency=$(squeue --noheader --format %i --name paper-run-indtest) \
  --wrap="source ~/.bashrc;cd /home/nrlab/wang04/ulyses/; git checkout logistic_reg;git pull;git add .;git commit -m 'paper run';git push"

