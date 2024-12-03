library(tidyverse)
library(openxlsx)

# set wd to /home/nrlab/wang04/ulyses/
setwd("/home/nrlab/wang04/ulyses/")
git_branch_name <- "logistic_reg"
git_command <- paste0("git checkout ", git_branch_name)
system(git_command)

out_excel_file <- "/home/nrlab/wang04/ulyses/6_statistics_source_data/Statistics Source Data.xlsx"

file01 <- "STATS_lr"
file01_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/STATS_lr_performance_of_all_feat_comb.csv"
file01_description <- "detailed stats of logistic regression models"

file02 <- "RAW_lr"
file02_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/lr_ranking_plot_main_fig.csv"
file02_description <- "raw performance of logistic regression model"

file03 <- "lr_indtest_score"
file03_on_server <- "/home/nrlab/wang04/ulyses/2_indtest_plots/lr_indtest_scores.csv"
file03_description <- "raw LR predicted scores in the test dataset"


file1 <- "STATS_xgb_x1-x6"
file1_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/xgboost_ranking_plot_main_fig.stats.csv"
file1_description <- "detailed stats related to Fig.5f and h; Fig.6a; Extended Data Fig. 3a; Extended Data Fig. 8-10a"

file2 <- "RAW_xgb_x1-x6"
file2_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/xgboost_ranking_plot_main_fig.csv"
file2_description <- "raw individual performance of x1 to x6 models (various performance metrics)"

file3 <- "STATS_xgb_all_feat"
file3_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/STATS_xgb_performance_of_all_feat_comb.csv"
file3_description <- "stats of numbers related to upset plot in Extended Data Fig. 3b-e"


file4 <- "RAW_xgb_all_feat"
file4_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/RAW_xgb_performance_of_all_feat_comb.csv"
file4_description <- "related to upset plot in Extended Data Fig. 3b-e"

file5 <- "xgb_indtest_score"
file5_on_server <- "/home/nrlab/wang04/ulyses/2_indtest_plots/xgboost_indtest_scores.csv"
file5_description <- "raw scores related to Fig. 6d, Extended Data Fig. 8d, Extended Data Fig. 9b and Extended Data Fig. 10d"

file6 <- "STATS_cnn_c1-c5"
file6_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/resnet18/C4/STATS_cnn_performance_c1-c5.csv"
file6_description <- "detailed stats related to Fig.5g-i; Fig. 6a; Extended Data Fig. 4; Extended Data Fig. 8-10a"

file7 <- "RAW_cnn_c1-c5"
file7_on_server <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/resnet18/C4/cnn_ranking_plot_main_fig.csv"
file7_description <- "raw individual performance of c1 to c6 models (various performance metrics)"

file8 <- "cnn_indtest_score"
file8_on_server <- "/home/nrlab/wang04/ulyses/2_indtest_plots/resnet18/C4/cnn_indtest_scores.csv"
file8_description <- "raw cnn predicted scores related to Fig. 6e, Extended Data Fig. 8e, Extended Data Fig. 9c and Extended Data Fig. 10e"

# command to run re-generate the files
file0_command <- "source ~/.bashrc;conda activate R4_3; Rscript /home/nrlab/wang04/ulyses/models_logistic_regression/lr_binary_performance_ranking.R;Rscript /home/nrlab/wang04/ulyses/models_logistic_regression/lr_binary_upset_plot.R"
file1_command <- "source ~/.bashrc;conda activate R4_3;Rscript /home/nrlab/wang04/ulyses/xgboost/xgboost_binary_performance_ranking.R"
file2_command <- file1_command
file5_command <- file1_command

file3_command <- "source ~/.bashrc;conda activate R4_3;Rscript /home/nrlab/wang04/ulyses/xgboost/xgboost_binary_upset_plot.R"
file4_command <- file3_command

file6_command <- "source ~/.bashrc;conda activate R4_3;Rscript /home/nrlab/wang04/ulyses/models/cnn_binary_performance_ranking_plot.R"
file7_command <- file6_command
file8_command <- file6_command

# sbatch commands

file0_submit_command <- paste0(
  'sbatch --job-name=stats_LR --output=stats_lr.out --error=stats_lr.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file0_command, '"'
)

file1_submit_command <- paste0(
  'sbatch --job-name=stats_xgb_performance_x1-x6 --output=stats_xgb_performance_x1-x6.out --error=stats_xgb_performance_x1-x6.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file1_command, '"'
)
file2_submit_command <- paste0(
  'sbatch --job-name=raw_xgb_performance_x1-x6 --output=raw_xgb_performance_x1-x6.out --error=raw_xgb_performance_x1-x6.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file2_command, '"'
)
file3_submit_command <- paste0(
  'sbatch --job-name=stats_xgb_performance_of_all_feat_comb --output=stats_xgb_performance_of_all_feat_comb.out --error=stats_xgb_performance_of_all_feat_comb.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file3_command, '"'
)
file4_submit_command <- paste0(
  'sbatch --job-name=raw_xgb_performance_of_all_feat_comb --output=raw_xgb_performance_of_all_feat_comb.out --error=raw_xgb_performance_of_all_feat_comb.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file4_command, '"'
)
file5_submit_command <- paste0(
  'sbatch --job-name=xgb_indtest_scores --output=xgb_indtest_scores.out --error=xgb_indtest_scores.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file5_command, '"'
)
file6_submit_command <- paste0(
  'sbatch --job-name=stats_cnn_performance_c1-c5 --output=stats_cnn_performance_c1-c5.out --error=stats_cnn_performance_c1-c5.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file6_command, '"'
)
file7_submit_command <- paste0(
  'sbatch --job-name=raw_cnn_performance_c1-c5 --output=raw_cnn_performance_c1-c5.out --error=raw_cnn_performance_c1-c5.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file7_command, '"'
)
file8_submit_command <- paste0(
  'sbatch --job-name=cnn_indtest_score --output=cnn_indtest_score.out --error=cnn_indtest_score.err --time=0-1 --mem=8G --partition=epyc --wrap="',
  file8_command, '"'
)


# run the commands
system(file0_command)
system(file1_command)
# system(file2_command)
system(file3_command)
# system(file4_command)
# system(file5_command)
system(file6_command)
# system(file7_command)
# system(file8_command)

# read in the data
file_toc_data <- tibble(
  file = c(
    file01,
    file02,
    file03,
    file1,
    file2,
    file3,
    file4,
    file5,
    file6,
    file7,
    file8
  ),
  description = c(
    file01_description,
    file02_description,
    file03_description,
    file1_description,
    file2_description,
    file3_description,
    file4_description,
    file5_description,
    file6_description,
    file7_description,
    file8_description
  )
)

file01_data <- read_csv(file01_on_server)
file02_data <- read_csv(file02_on_server)
file03_data <- read_csv(file03_on_server)
file1_data <- read_csv(file1_on_server)
file2_data <- read_csv(file2_on_server)
file3_data <- read_csv(file3_on_server)
file4_data <- read_csv(file4_on_server)
file5_data <- read_csv(file5_on_server)
file6_data <- read_csv(file6_on_server)
file7_data <- read_csv(file7_on_server)
file8_data <- read_csv(file8_on_server)

# write to excel as different sheets

# Create a list of tibbles
tibbles <- list(
  Index = file_toc_data,
  file01 = file01_data,
  file02 = file02_data,
  file03 = file03_data,
  file1 = file1_data,
  file2 = file2_data,
  file3 = file3_data,
  file4 = file4_data,
  file5 = file5_data,
  file6 = file6_data,
  file7 = file7_data,
  file8 = file8_data
)

names(tibbles) <- c(
  "Table of Contents",
  file01,
  file02,
  file03,
  file1,
  file2,
  file3,
  file4,
  file5,
  file6,
  file7,
  file8
)

# Write the list of tibbles to an Excel file
write.xlsx(x = tibbles, file = out_excel_file)
message("Data written to ", out_excel_file)
