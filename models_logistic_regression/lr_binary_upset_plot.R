library(ggupset)
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(argparse)
library(kableExtra)
library(gridExtra)
library(grid)
library(nord)
library(ggpubr)



# add parameters
parser <- argparse::ArgumentParser()
# add tf_label argument, default to 'all'
parser$add_argument("--tf_label", type = "character", help = "tf_label", default = "0-0.03_pairwise")
# add wd argument, default as current directory
parser$add_argument("--wd", type = "character", help = "wd", default = "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/lr_model/lr10")
# add plot_width argument, default to 80
parser$add_argument("--plot_width", type = "integer", help = "plot_width", default = 85)
# add plot_height argument, default to 170
parser$add_argument("--plot_height", type = "integer", help = "plot_height", default = 120)
# add plot_unit argument, default to 'mm'
parser$add_argument("--plot_unit", type = "character", help = "plot_unit", default = "mm")
# add plot_dpi argument, default to 300
parser$add_argument("--plot_dpi", type = "integer", help = "plot_dpi", default = 300)
# add n_repeat argument, default to 50
parser$add_argument("--n_repeat", type = "integer", help = "n_repeat", default = 50)
# add all_feat_name, default as "cnv_sd_length_motif"
parser$add_argument("--all_feat_name", type = "character", help = "all_feat_name", default = "cnv_sd_length_ctRatio_slRatio")
# add all_feat_name_new, default as "Length+SD+CNV+C/T+S/L"
parser$add_argument("--all_feat_name_new", type = "character", help = "all_feat_name_new", default = "All features")
parser$add_argument("--manual_colors", type = "character", help = "manual_colors", default = c("#8fbcbb", "#d088c0", "#e0869a", "#eed093", "#8f8cd4"))

# parse args
args <- parser$parse_args()
tf_label <- args$tf_label
wd <- args$wd
plot_width <- args$plot_width
plot_height <- args$plot_height
plot_unit <- args$plot_unit
plot_dpi <- args$plot_dpi
n_repeat <- args$n_repeat
all_feat_name <- args$all_feat_name
all_feat_name_new <- args$all_feat_name_new
manual_colors <- args$manual_colors

nested_cv_plot_path <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots"

source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

# based on tf_label, set the y axis label
tf_range <- if (tf_label == "0.01-0.03_pairwise") {
  "(0.01, 0.03]"
} else if (tf_label == "0.03-0.1_pairwise") {
  "(0.03, 0.1]"
} else if (tf_label == "0.1-1_pairwise") {
  "(0.1, 1]"
} else if (tf_label == "0.2-1_pairwise") {
  "(0.2, 1]"
} else if (tf_label == "0-0.03_pairwise") {
  "(0, 0.03]"
} else if (tf_label == "all") {
  "all"
}




plot_file <- file.path(wd, tf_label, paste("final_boxplot_tf_", tf_label, ".pdf", sep = ""))
###############################################################################
# plot the performance
###############################################################################
fl <- list.files(
  path = wd,
  "repeat_10_fold_5_metrics.csv", recursive = TRUE, full.names = TRUE
)

fl_read <- bettermc::mclapply(fl, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

fl_read <- lapply(fl, read_csv, show_col_types = FALSE)
names(fl_read) <- fl
fl_read_bind <- bind_rows(fl_read, .id = "id")

fl_read_bind <- bind_rows(fl_read, .id = "id") %>%
  mutate(ichorcna_strat = case_when(
    stringr::str_detect(id, "0-0\\.03") ~ "[0, 0.03]",
    stringr::str_detect(id, "0\\.03-0\\.1") ~ "(0.03, 0.1]",
    stringr::str_detect(id, "0\\.1-1") ~ "(0.1, 1]",
    stringr::str_detect(id, "all") ~ "all",
    TRUE ~ "unknown"
  ))
# set order of ichorcna_strat
fl_read_bind$ichorcna_strat <- factor(fl_read_bind$ichorcna_strat,
  levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all")
)

# make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
fl_read_bind <- fl_read_bind %>%
  mutate(feat = "TF")



# plot the feature performance of 'cnv_sd_length_motif'
ans2 <- fl_read_bind %>%
  # filter(feat %in% c(all_feat_name, "length", "sd", "cnv", "slRatio", "ctRatio")) %>%
  # filter out ichorcna_strat == "unknown"
  filter(ichorcna_strat != "unknown") %>%
  # remove where the value of test_threshold_99spe is Inf
  # filter(!is.infinite(test_threshold_99spe)) %>%
  # remove where the value of test_threshold_98spe is Inf
  # filter(!is.infinite(test_threshold_98spe)) %>%
  # remove where the value of test_threshold_95spe is Inf
  # filter(!is.infinite(test_threshold_95spe)) %>%
  # pivot_longer to make .metric as a column, all cols start with 'test_' to ".metric", values to ".estimate"
  pivot_longer(cols = starts_with("test_"), names_to = ".metric", values_to = ".estimate") %>%
  # calculate the mean .metric of each feature of each ichorcna_strat
  # group_by(ichorcna_strat, feat, .metric) %>%
  # summarize(mean = mean(.estimate), median = median(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  # change 'roc_auc' to 'AUROC'
  mutate(.metric = ifelse(.metric == "roc_auc", "AUROC", .metric)) %>%
  # change 'accuracy' to 'Accuracy'
  mutate(.metric = ifelse(.metric == "accuracy", "Accuracy", .metric))
# # change 'slRatio' to 'S/L'
# mutate(feat = ifelse(feat == "slRatio", "S/L", feat)) %>%
# # change 'ctRatio' to 'C/T'
# mutate(feat = ifelse(feat == "ctRatio", "C/T", feat)) %>%
# # change 'cnv' to 'CNV'
# mutate(feat = ifelse(feat == "cnv", "CNV", feat)) %>%
# # change 'sd' to 'SD'
# mutate(feat = ifelse(feat == "sd", "SD", feat)) %>%
# # change 'length' to 'Length'
# mutate(feat = ifelse(feat == "length", "Length", feat)) %>%
# # change 'all_feat_name' to 'Length+SD+CNV+C/T+S/L'
# mutate(feat = ifelse(feat == all_feat_name, all_feat_name_new, feat))

# set the order of feat
manual_colors2 <- c("grey3", manual_colors[1:length(unique(ans2$feat))])



ans2 <- ans2 %>%
  filter(!str_detect(.metric, "threshold"))

ans2_without_all <- ans2 %>%
  filter(ichorcna_strat != "all") %>%
  filter(!str_detect(.metric, "threshold"))

# fl_read_bind2 <- ans2 %>%
#   # add Length column
#   mutate("Length" = ifelse(str_detect(feat, "length"), TRUE, FALSE)) %>%
#   # add SD column
#   mutate("SD" = ifelse(str_detect(feat, "sd"), TRUE, FALSE)) %>%
#   # add CNV column
#   mutate("CNV" = ifelse(str_detect(feat, "cnv"), TRUE, FALSE)) %>%
#   # add S/L Ratio column
#   mutate("S/L" = ifelse(str_detect(feat, "slRatio"), TRUE, FALSE)) %>%
#   # add C/T Ratio column
#   mutate("C/T" = ifelse(str_detect(feat, "ctRatio"), TRUE, FALSE))


# tmp_func <- function(Length, SD, CNV, `S/L`, `C/T`) {
#   c(if (Length) "Length", if (SD) "SD", if (CNV) "CNV", if (`S/L`) "S/L", if (`C/T`) "C/T")
# }

# fl_read_bind3 <- fl_read_bind2 %>%
#   # add Length column
#   mutate("Label" = pmap(list(Length, SD, CNV, `S/L`, `C/T`), tmp_func))

# plot ichorcna_strat == [0, 0.03], .metric == "test_sen_99spe"

# write a df to expand the all combinations of 'ichorcna_strat' and '.metric'
# df <- expand.grid(ichorcna_strat_param = unique(ans2$ichorcna_strat), .metric_param = unique(ans2$.metric))
# # make a new col called 'ichorcna_strat_label', substitute all blanks with '_'
# df$ichorcna_strat_label <- gsub(" ", "_", df$ichorcna_strat_param)
# df <- df %>%
#   mutate(plot_file = paste0(nested_cv_plot_path, "/lr_upsetplot_tf_", ichorcna_strat_label, "_", .metric_param, ".pdf")) %>%
#   select(-ichorcna_strat_label)


# map the function to the df
# pmap_df(df, upset_plot2, input_df = fl_read_bind3)

# save fl_read_bind3 as rds file
# saveRDS(fl_read_bind3, file = paste0(nested_cv_plot_path, "/final_upsetplot_input_data", ".rds"))

# save fl_read_bind3 as csv file
write_csv(ans2, paste0(nested_cv_plot_path, "/RAW_lr_performance_of_all_feat_comb", ".csv"))


fl_read_bind3_stats <- ans2 %>%
  group_by(ichorcna_strat, feat, .metric) %>%
  summarize(
    mean = Hmisc::smean.cl.boot(.estimate)["Mean"],
    median = median(.estimate),
    sd = Hmisc::smean.sd(.estimate)["SD"],
    ci_95_lower = Hmisc::smean.cl.boot(.estimate)["Lower"],
    ci_95_upper = Hmisc::smean.cl.boot(.estimate)["Upper"],
    sem_lower = mean - sd / sqrt(n()),
    sem_upper = mean + sd / sqrt(n()),
    .groups = "drop"
  )
# summarize(
#   mean = mean(.estimate),
#   median = median(.estimate),
#   sd = sd(.estimate),
#   # lower 95CI
#   ci_95_lower = mean - 1.96 * sd / sqrt(n()),
#   # upper 95CI
#   ci_95_upper = mean + 1.96 * sd / sqrt(n()),
#   # standard error of mean lower
#   sem_lower = mean - sd / sqrt(n()),
#   # standard error of mean upper
#   sem_upper = mean + sd / sqrt(n()),
#   .groups = "drop"
# )


write_csv(fl_read_bind3_stats, paste0(nested_cv_plot_path, "/STATS_lr_performance_of_all_feat_comb", ".csv"))
message("Saved the stats of the performance of all features in ", nested_cv_plot_path, "/STATS_lr_performance_of_all_feat_comb.csv")
