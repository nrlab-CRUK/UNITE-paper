library(ggupset)
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(argparse)
library(kableExtra)
library(gridExtra)
library(grid)
library(nord)
# library(gghalves) and install it if not exists
if (!requireNamespace("gghalves", quietly = TRUE)) {
  install.packages("gghalves")
}

library(gghalves)
###############################################################################
# add parameters
###############################################################################

# add parameters
parser <- argparse::ArgumentParser()
# add wd argument, default as current directory
parser$add_argument("--wd", type = "character", help = "wd", default = "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/lr_model/lr10")
# add all_feat_name, default as "cnv_sd_length_motif"
parser$add_argument("--all_feat_name", type = "character", help = "all_feat_name", default = "cnv_sd_length_ctRatio_slRatio")
# add all_feat_name_new, default as "Length+SD+CNV+C/T+S/L"
parser$add_argument("--all_feat_name_new", type = "character", help = "all_feat_name_new", default = "All")
# add manual_colors, default as c("#8fbcbb", "#d088c0", "#e0869a", "#eed197", "#8f8cd4")
parser$add_argument("--manual_colors", type = "character", help = "manual_colors", default = c("#8fbcbb", "#d088c0", "#e0869a", "#eed093", "#8f8cd4"))
# add plot_width, default as 200
parser$add_argument("--plot_width", type = "numeric", help = "plot_width", default = 120)
# add plot_height, default as 100
parser$add_argument("--plot_height", type = "numeric", help = "plot_height", default = 50)
# add plot_unit, default as "mm"
parser$add_argument("--plot_unit", type = "character", help = "plot_unit", default = "mm")
# add labs_x, default as "ichorCNA TF Strata"
parser$add_argument("--labs_x", type = "character", help = "labs_x", default = NULL)
# add labs_y, default as "Mean Performance"
parser$add_argument("--labs_y", type = "character", help = "labs_y", default = "Mean performance and 95% CI")

# add all_hight_color, set default as "#b5dae6"
parser$add_argument("--all_hight_color", type = "character", help = "all_hight_color", default = "#c1dde7")
# add metrics_shown_in_main_fig
parser$add_argument("--metrics_shown_in_main_fig", type = "character", help = "metrics_shown_in_main_fig", default = c("AUROC"))

# add a param called "what_to_show_in_appendix", default as c("AUROC", "AUPRC", "Sensitivity", "Specificity", "F1", "Accuracy", "PPV")
parser$add_argument("--what_to_show_in_appendix", type = "character", help = "what_to_show_in_appendix", default = c("AUROC", "AUPRC", "Sensitivity", "Specificity", "F1", "Accuracy", "PPV"))


# parse args
args <- parser$parse_args()
wd <- args$wd
all_feat_name <- args$all_feat_name
all_feat_name_new <- args$all_feat_name_new
manual_colors <- args$manual_colors
plot_width <- args$plot_width
plot_height <- args$plot_height
plot_unit <- args$plot_unit
labs_x <- args$labs_x
labs_y <- args$labs_y
all_hight_color <- args$all_hight_color
metrics_shown_in_main_fig <- args$metrics_shown_in_main_fig
what_to_show_in_appendix <- args$what_to_show_in_appendix

nested_cv_plot_path <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots"
indtest_plot_path <- "/home/nrlab/wang04/ulyses/2_indtest_plots"
plot_file <- file.path(nested_cv_plot_path, "lr_ranking_plot_main_fig.pdf")
plot_file_csv <- file.path(nested_cv_plot_path, "lr_ranking_plot_main_fig.csv")
plot_file_indtest <- file.path(indtest_plot_path, "lr_ranking_plot_indtest_main_fig.pdf")
indtest_csv <- file.path(indtest_plot_path, "lr_indtest_scores.csv")
plot_file_plus_all_feat <- file.path(nested_cv_plot_path, "lr_ranking_plot_all_feat_appendix.pdf")
plot_file_plus_all_feat_indtest <- file.path(indtest_plot_path, "lr_ranking_plot_all_feat_indtest_appendix.pdf")




source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

###############################################################################
# plot the performance
###############################################################################

fl <- list.files(
  path = wd,
  "repeat_10_fold_5_metrics.csv", recursive = TRUE, full.names = TRUE
)

fl_read <- bettermc::mclapply(fl, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_read) <- fl

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
  pivot_longer(cols = starts_with("test_"), names_to = ".metric", values_to = ".estimate")


# ans2 <- ans2 %>%
# calculate the mean .metric of each feature of each ichorcna_strat
# group_by(ichorcna_strat, feat, .metric) %>%
# summarize(mean = mean(.estimate), median = median(.estimate), sd = sd(.estimate), .groups = "drop") %>%
# change 'slRatio' to 'S/L'
# mutate(feat = ifelse(feat == "slRatio", "S/L", feat)) %>%
# change 'ctRatio' to 'C/T'
# mutate(feat = ifelse(feat == "ctRatio", "C/T", feat)) %>%
# change 'cnv' to 'CNV'
# mutate(feat = ifelse(feat == "cnv", "CNV", feat)) %>%
# change 'sd' to 'SD'
# mutate(feat = ifelse(feat == "sd", "SD", feat)) %>%
# change 'length' to 'Length'
# mutate(feat = ifelse(feat == "length", "Length", feat)) %>%
# change 'all_feat_name' to 'Length+SD+CNV+C/T+S/L'
# mutate(feat = ifelse(feat == all_feat_name, all_feat_name_new, feat))

# set the order of feat
feat_lv <- c("TF")
ans2$feat <- factor(ans2$feat, levels = feat_lv)
lr_model_label <- c("LR")
mainplot_legend_labels <- paste(lr_model_label, feat_lv, sep = "-")
names(mainplot_legend_labels) <- feat_lv
manual_colors2 <- c("grey3", manual_colors[1:length(unique(ans2$feat))])



ans2 <- ans2 %>%
  filter(!str_detect(.metric, "test_threshold")) %>%
  # filter(!(.metric %in% c("test_sensitivity", "test_specificity", "test_precision", "test_recall", "test_f1", "test_acc", "test_ppv")))
  filter(!(.metric %in% c("test_precision", "test_recall")))

# set the factor levels of ans2$.metric
ans2$.metric <- factor(ans2$.metric, levels = c(
  "test_auroc",
  "test_auprc",
  "test_sensitivity",
  "test_specificity",
  "test_f1",
  "test_acc",
  "test_ppv",
  "test_sen_99spe", "test_sen_98spe", "test_sen_95spe",
  "test_acc_99spe", "test_acc_98spe", "test_acc_95spe",
  "test_f1_99spe", "test_f1_98spe", "test_f1_95spe",
  "test_ppv_threshold_99spe", "test_ppv_threshold_98spe", "test_ppv_threshold_95spe"
))

# rename .metric
ans2$.metric <- recode(ans2$.metric,
  "test_auroc" = "AUROC",
  "test_auprc" = "AUPRC",
  "test_sensitivity" = "Sensitivity",
  "test_specificity" = "Specificity",
  "test_f1" = "F1",
  "test_acc" = "Accuracy",
  "test_ppv" = "PPV",
  "test_sen_99spe" = "Sensitivity at 99% Specificity",
  "test_sen_98spe" = "Sensitivity at 98% Specificity",
  "test_sen_95spe" = "Sensitivity at 95% Specificity",
  "test_acc_99spe" = "Accuracy at 99% Specificity",
  "test_acc_98spe" = "Accuracy at 98% Specificity",
  "test_acc_95spe" = "Accuracy at 95% Specificity",
  "test_f1_99spe" = "F1 at 99% Specificity",
  "test_f1_98spe" = "F1 at 98% Specificity",
  "test_f1_95spe" = "F1 at 95% Specificity",
  "test_ppv_threshold_99spe" = "PPV at 99% Specificity",
  "test_ppv_threshold_98spe" = "PPV at 98% Specificity",
  "test_ppv_threshold_95spe" = "PPV at 95% Specificity"
)

# remove .metric contain the text "PPV"
ans2 <- ans2 %>%
  filter(!str_detect(.metric, "PPV")) %>%
  filter(!str_detect(.metric, "98%"))
pd <- position_dodge(0.1)
###############################################################################
# appendix plot
###############################################################################

ans2_appendix <- dplyr::filter(ans2, .metric %in% what_to_show_in_appendix) %>%
  mutate(.metric = factor(.metric, levels = what_to_show_in_appendix))
# do a plot, x axis is ichorcna_strat, y axis is mean .metric, color is feat, facet by .metric
p_appendix <- ggplot(data = ans2_appendix, aes(x = ichorcna_strat, y = .estimate, color = feat, fill = feat, group = feat)) +
  # Add background color to 'all' group
  geom_rect(
    data = subset(ans2_appendix, ichorcna_strat == "all"),
    aes(
      xmin = as.numeric(ichorcna_strat) - 0.5,
      xmax = as.numeric(ichorcna_strat) + 0.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = all_hight_color,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  # make dot size to 0.8 and add outer line
  stat_summary(fun = mean, geom = "line", data = ans2_appendix %>% filter(ichorcna_strat != "all")) +
  stat_summary(fun = mean, geom = "point", shape = 21, color = "grey2", size = 1) +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "errorbar",
    # color = "black",
    position = pd,
    linewidth = 0.1
  ) +
  # geom_line(data = ans2_without_all, aes(x = ichorcna_strat, y = mean, color = feat, group = feat)) +
  # geom_point(size = 0.8) +
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  facet_wrap(~.metric, scales = "free_y", drop = TRUE) +
  # set the color of feat to high contrast but elegant colors
  scale_color_manual(values = manual_colors2) +
  scale_fill_manual(values = manual_colors2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # add legend
  theme(strip.background = element_rect(fill = "#e9e6e6", color = "#e9e6e6")) +
  labs(x = labs_x, y = labs_y) +
  # remove x-axis coordiante extension
  scale_x_discrete(expand = c(0.05, 0)) +
  theme(panel.spacing = unit(0.1, "cm")) +
  theme(strip.text = element_text(size = 5)) +
  theme(axis.text = element_text(size = 5)) +
  # remove legend title
  theme(axis.title = element_text(size = 6)) +
  theme(legend.position = "top") +
  theme(legend.text = element_text(size = 5), legend.title = element_blank()) +
  theme(legend.key = element_rect(fill = "transparent", color = "transparent")) +
  theme(legend.background = element_rect(fill = "transparent", color = "transparent")) +
  theme(legend.box.background = element_rect(fill = "transparent", color = "transparent"))

p_appendix_no_legend <- p_appendix + theme(legend.position = "none")

ggsave(plot_file_plus_all_feat,
  p_appendix,
  create.dir = TRUE,
  width = 180, height = 100, unit = plot_unit
)
message(plot_file_plus_all_feat)

plot_file_plus_all_feat_no_legend <- gsub(".pdf", "_no_legend.pdf", plot_file_plus_all_feat)
ggsave(plot_file_plus_all_feat_no_legend,
  p_appendix_no_legend,
  create.dir = TRUE,
  width = 180, height = 90, unit = plot_unit
)

message(plot_file_plus_all_feat)


###############################################################################
# main figure
###############################################################################
ans2_main <- ans2 %>%
  filter(.metric %in% metrics_shown_in_main_fig) %>%
  mutate(.metric = factor(.metric, levels = metrics_shown_in_main_fig))


p_main <- ggplot(data = ans2_main, aes(x = ichorcna_strat, y = .estimate, color = feat, fill = feat, group = feat)) +
  # Add background color to 'all' group
  geom_rect(
    data = subset(ans2_main, ichorcna_strat == "all"),
    aes(
      xmin = as.numeric(ichorcna_strat) - 0.5,
      xmax = as.numeric(ichorcna_strat) + 0.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = all_hight_color,
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  # make dot size to 0.8 and add outer line
  # geom_point(data = ans2, size = 0.99, shape = 21, fill = "transparent", color = "grey2") +
  # geom_point(data = ans2, size = 0.99, shape = 21, fill = feat, color = "grey2") +
  stat_summary(fun = mean, geom = "line", data = ans2_main %>% filter(ichorcna_strat != "all"), linewidth = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 21, color = "grey2", size = 2) +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "errorbar",
    # color = "black",
    position = pd,
    linewidth = 0.3
  ) +
  ggpubr::stat_compare_means(
    method = "kruskal", label = "p.signif", label.y = c(0.91, 0.96, 0.993, 0.93), size = 3, color = "#756a6a",
    paired = FALSE
  ) +
  # geom_line(data = ans2_without_all, aes(x = ichorcna_strat, y = mean, color = feat, group = feat)) +
  # geom_point(size = 0.8) +
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  facet_wrap(~.metric, scales = "free_y") +
  # set the color of feat to high contrast but elegant colors
  scale_color_manual(values = manual_colors2, labels = mainplot_legend_labels) +
  scale_fill_manual(values = manual_colors2, labels = mainplot_legend_labels) +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  # add legend
  theme(strip.background = element_rect(fill = "#e9e6e6", color = "#e9e6e6")) +
  labs(x = labs_x, y = "Mean performance and 95% CI") +
  # remove x-axis coordiante extension
  scale_x_discrete(expand = c(0.05, 0)) +
  theme(panel.spacing = unit(0.1, "cm")) +
  theme(strip.text = element_text(size = 6)) +
  theme(axis.text = element_text(size = 6)) +
  # remove legend title
  theme(axis.title = element_text(size = 6)) +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size = 6), legend.title = element_blank()) +
  theme(legend.key = element_rect(fill = "transparent", color = "transparent")) +
  theme(legend.background = element_rect(fill = "transparent", color = "transparent")) +
  theme(
    legend.box.background = element_rect(fill = "transparent", color = "transparent"),
    legend.box.spacing = unit(0, "pt"), # The spacing between the plotting area and the legend box (unit)
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0)
  )



ggsave(plot_file,
  p_main,
  width = 63,
  height = 60,
  unit = plot_unit,
  create.dir = TRUE
)

message(plot_file)

# save the plot to rds file
saveRDS(p_main, gsub(".pdf", ".rds", plot_file))

# hide legend
p_main_no_legend <- p_main + theme(legend.position = "none")
plot_file_no_legend <- gsub(".pdf", "_no_legend.pdf", plot_file)
ggsave(plot_file_no_legend,
  p_main_no_legend,
  width = 50,
  height = 60,
  unit = plot_unit,
  create.dir = TRUE
)

message(plot_file_no_legend)
# save the raw data to a csv
ans2 %>%
  write_csv(plot_file_csv)

message("Save the raw data to: ", plot_file_csv)


# for each ichorcna_strat, feat, .metric, calculate the mean, median, sd, 95% CI
ans2_appendix_stats <- ans2_appendix %>%
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


# save the stats to a csv
ans2_appendix_stats %>%
  write_csv(gsub(".csv", ".stats.csv", plot_file_csv))


###############################################################################
# plot dots of the indtest scores
###############################################################################


# get indtest scores

fl_indtest <- list.files(
  path = wd,
  "repeat_0_fold_0_metrics.csv", recursive = TRUE, full.names = TRUE
)

fl_indtest_read <- bettermc::mclapply(fl_indtest, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_indtest_read) <- fl_indtest

fl_indtest_read_bind <- bind_rows(fl_indtest_read, .id = "id") %>%
  mutate(ichorcna_strat = case_when(
    stringr::str_detect(id, "0-0\\.03") ~ "[0, 0.03]",
    stringr::str_detect(id, "0\\.03-0\\.1") ~ "(0.03, 0.1]",
    stringr::str_detect(id, "0\\.1-1") ~ "(0.1, 1]",
    stringr::str_detect(id, "all") ~ "all",
    TRUE ~ "unknown"
  ))
# set order of ichorcna_strat
fl_indtest_read_bind$ichorcna_strat <- factor(fl_indtest_read_bind$ichorcna_strat,
  levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all")
)

# make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
fl_indtest_read_bind <- fl_indtest_read_bind %>%
  mutate(feat = "TF")




# plot the feature performance of 'cnv_sd_length_motif'
fl_indtest_read_bind <- fl_indtest_read_bind %>%
  # filter(feat %in% c(all_feat_name, "length", "sd", "cnv", "slRatio", "ctRatio")) %>%
  # filter out ichorcna_strat == "unknown"
  filter(ichorcna_strat != "unknown") %>%
  # test_sen_98spe is 1 when test_auroc is 1
  mutate(test_sen_98spe = ifelse(test_auroc == 1, 1, test_sen_98spe)) %>%
  # pivot_longer to make .metric as a column, all cols start with 'test_' to ".metric", values to ".estimate"
  pivot_longer(cols = starts_with("test_"), names_to = ".metric", values_to = ".estimate")

fl_indtest_read_bind_all_feature_model <- fl_indtest_read_bind %>%
  filter(feat == "TF")

test_99spec_thresholds <- fl_indtest_read_bind_all_feature_model %>%
  filter(.metric == "test_threshold_99spe") %>%
  select(ichorcna_strat, .metric, .estimate) %>%
  distinct()

test_95spec_thresholds <- fl_indtest_read_bind_all_feature_model %>%
  filter(.metric == "test_threshold_95spe") %>%
  select(ichorcna_strat, .metric, .estimate) %>%
  distinct()

thresholds <- fl_indtest_read_bind_all_feature_model %>%
  select(ichorcna_strat, .metric, .estimate) %>%
  distinct()

###############################################################################

meta_data <- read_csv(meta_csv_latest)

fl <- list.files(
  path = wd,
  "repeat_0_fold_0_y_test_y_pred_z_test.csv", recursive = TRUE, full.names = TRUE
)

fl_read <- bettermc::mclapply(fl, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_read) <- fl

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


# save fl_read_bind to a csv
fl_read_bind %>%
  write_csv(file.path(wd, "fl_read_bind_indtest.csv"))


# only keep the all_feature model, by str_detect(id, "cnv_sd_length_ctRatio_slRatio")
fl_read_bind <- fl_read_bind %>%
  filter(str_detect(feat, "TF"))


# add meta_data$cohort to ans2, based on the bam_id
meta_data_selected <- meta_data %>%
  select(bam_id, cohort)

ans2 <- left_join(fl_read_bind, meta_data_selected, by = c("bam_id_test" = "bam_id"))

# y_test as factor
ans2$y_test <- factor(ans2$y_test, levels = c(0, 1))

# recode the y_test to "Healthy" and "Cancer"
ans2$y_test <- recode(ans2$y_test, "0" = "Healthy", "1" = "Cancer")

# make the cohort as factor
ans2$cohort <- factor(ans2$cohort)

# set the level of cohort, make "Healthy" the first level using relevel
ans2$cohort <- relevel(ans2$cohort, ref = "Healthy")

# report the number of samples in each ichorcna_strat and cohort
ans2_summary <- ans2 %>%
  group_by(ichorcna_strat, cohort) %>%
  summarize(n = n())

df_to_save <- ans2
df <- expand.grid(ichorcna_strat_param = unique(ans2$ichorcna_strat))
# make a new col called 'ichorcna_strat_label', substitute all blanks with '_'
df$ichorcna_strat_label <- gsub(" ", "_", df$ichorcna_strat_param)
df <- df %>%
  mutate(plot_file = paste0(indtest_plot_path, "/lr_pred_score_plot_indtest_tf_", ichorcna_strat_label, ".pdf")) %>%
  select(-ichorcna_strat_label)

score_plot <- function(ichorcna_strat_param, plot_file, thresholds, ans2, point_shape = 21) {
  # plot
  target_ichorcna_strat <- ichorcna_strat_param
  target_99threshold <- thresholds %>%
    filter(.metric == "test_threshold_99spe") %>%
    filter(ichorcna_strat == target_ichorcna_strat) %>%
    pull(.estimate)

  target_95threshold <- thresholds %>%
    filter(.metric == "test_threshold_95spe") %>%
    filter(ichorcna_strat == target_ichorcna_strat) %>%
    pull(.estimate)

  sen_99spec <- thresholds %>%
    filter(.metric == "test_sen_99spe") %>%
    filter(ichorcna_strat == target_ichorcna_strat) %>%
    pull(.estimate)

  sen_99_spec_label <- sen_99spec %>%
    round(., 4) %>%
    # convert to percentage and append %
    scales::percent(., accuracy = 0.01) %>%
    paste0("Sen = ", .)
  sen_95spec <- thresholds %>%
    filter(.metric == "test_sen_95spe") %>%
    filter(ichorcna_strat == target_ichorcna_strat) %>%
    pull(.estimate)
  sen_95_spec_label <- sen_95spec %>%
    round(., 4) %>%
    # convert to percentage and append %
    scales::percent(., accuracy = 0.01) %>%
    paste0("Sen = ", .)


  y_label_for_99spec <- paste0("99% Overall Spe\n", "Threshold = ", round(target_99threshold, 4), "\n", sen_99_spec_label)
  y_label_for_95spec <- paste0("95% Overall Spe\n", "Threshold = ", round(target_95threshold, 4), "\n", sen_95_spec_label)

  # plot as dot plot, x axis is cohort, y axis is y_pred, color is y_test, facet by ichorcna_strat
  data <- ans2 %>%
    filter(ichorcna_strat == target_ichorcna_strat)

  # add a col called cohort_detected_by_99threshold, if y_pred > target_99threshold, then cohort_detected_by_99threshold is 1, else 0
  data <- data %>%
    mutate(cohort_detected_by_99threshold = ifelse(y_pred > target_99threshold, 1, 0))
  # add a col called cohort_detected_by_95threshold, if y_pred > target_99threshold, then cohort_detected_by_99threshold is 1, else 0
  data <- data %>%
    mutate(cohort_detected_by_95threshold = ifelse(y_pred > target_95threshold, 1, 0))

  # add a col called "cohort_sum_detected_by_99threshold", sum of cohort_detected_by_99threshold
  data <- data %>%
    group_by(cohort) %>%
    mutate(cohort_sum_detected_by_99threshold = sum(cohort_detected_by_99threshold)) %>%
    mutate(cohort_sum_detected_by_95threshold = sum(cohort_detected_by_95threshold)) %>%
    mutate(cohort_sum = n()) %>%
    ungroup()

  # add a col called cohort_sen_by_99threshold, cohort_sum_detected_by_99threshold / cohort_sum
  data <- data %>%
    mutate(cohort_sen_by_99spec = cohort_sum_detected_by_99threshold / cohort_sum) %>%
    mutate(cohort_sen_by_95spec = cohort_sum_detected_by_95threshold / cohort_sum)

  # add a col called cohort_sen_by_99spec_label, values to  paste0("sen99spec", "( ", cohort_sum_detected_by_99threshold, "/", cohort_sum, ",", cohort_sen_by_99spec, ")")
  data <- data %>%
    mutate(cohort_sen_by_99spec_label = case_when(
      cohort == "Healthy" ~ paste0("99% spe:", cohort_sum_detected_by_99threshold, "/", cohort_sum, ", ", round(cohort_sen_by_99spec, 4) %>% scales::percent(., accuracy = 0.01)),
      TRUE ~ paste0(cohort_sum_detected_by_99threshold, "/", cohort_sum, ", ", round(cohort_sen_by_99spec, 4) %>% scales::percent(., accuracy = 0.01))
    )) %>%
    mutate(cohort_sen_by_95spec_label = case_when(
      cohort == "Healthy" ~ paste0("95% spe:", cohort_sum_detected_by_95threshold, "/", cohort_sum, ", ", round(cohort_sen_by_95spec, 4) %>% scales::percent(., accuracy = 0.01)),
      TRUE ~ paste0(cohort_sum_detected_by_95threshold, "/", cohort_sum, ", ", round(cohort_sen_by_95spec, 4) %>% scales::percent(., accuracy = 0.01))
    ))

  # perfect case
  data <- data %>%
    mutate(cohort_sen_by_99spec_label_perfect = case_when(
      cohort == "Healthy" ~ paste0("99% spe:", "0", "/", cohort_sum, ", ", round(0, 4) %>% scales::percent(., accuracy = 0.01)),
      TRUE ~ paste0(cohort_sum, "/", cohort_sum, ", ", round(1, 4) %>% scales::percent(., accuracy = 0.01))
    )) %>%
    mutate(cohort_sen_by_95spec_label_perfect = case_when(
      cohort == "Healthy" ~ paste0("95% spe: ", "0", "/", cohort_sum, ", ", round(0, 4) %>% scales::percent(., accuracy = 0.01)),
      TRUE ~ paste0(cohort_sum, "/", cohort_sum, ", ", round(1, 4) %>% scales::percent(., accuracy = 0.01))
    ))

  # make cohort_label
  # not perfect case
  data <- data %>%
    mutate(cohort_label = paste0(cohort, "\n", cohort_sen_by_99spec_label, "\n", cohort_sen_by_95spec_label))
  healthy_label <- data %>%
    filter(str_detect(cohort_label, "Healthy")) %>%
    pull(cohort_label) %>%
    unique() %>%
    as.character()
  data$cohort_label <- as.factor(data$cohort_label)
  data$cohort_label <- relevel(data$cohort_label, ref = healthy_label)
  # perfect case
  data <- data %>%
    mutate(cohort_label_perfect = paste0(cohort, "\n", cohort_sen_by_99spec_label_perfect, "\n", cohort_sen_by_95spec_label_perfect))
  healthy_label_perfect <- data %>%
    filter(str_detect(cohort_label_perfect, "Healthy")) %>%
    pull(cohort_label_perfect) %>%
    unique() %>%
    as.character()
  data$cohort_label_perfect <- as.factor(data$cohort_label_perfect)
  data$cohort_label_perfect <- relevel(data$cohort_label_perfect, ref = healthy_label_perfect)


  # add a col called 'y_pred_99threshold' to data, if y_pred > target_99threshold, then 'y_pred_99threshold' is "Cancer", else "Healthy"
  data <- data %>%
    mutate(y_pred_99threshold = ifelse(y_pred > target_99threshold, "Cancer", "Healthy")) %>%
    mutate(y_pred_99threshold = factor(y_pred_99threshold, levels = c("Healthy", "Cancer")))

  # add a col called "dot_col", if y_pred > target_99threshold, then "tomato3", else if y_pred > target_95threshold & y_pred <= target_99threshold , then "orange", else "grey"
  data <- data %>%
    mutate(dot_col = case_when(
      y_pred > target_99threshold ~ "detected_by_99threshold",
      y_pred > target_95threshold & y_pred <= target_99threshold ~ "detected_by_95threshold",
      TRUE ~ "not_detected"
    ))

  dot_col_manual <- c("detected_by_99threshold" = "#251714FF", "detected_by_95threshold" = "#D9636CFF", "not_detected" = "lightgrey")
  # dot_col_manual <- c("detected_by_99threshold" = "#0d0887", "detected_by_95threshold" = "#cc4778", "not_detected" = "#f0f921")
  # dot_col_manual <- c("detected_by_99threshold" = "#0d0887", "detected_by_95threshold" = "#f0f921", "not_detected" = "lightgrey")
  if (target_ichorcna_strat != "(0.1, 1]") {
    p <- ggplot(data = data, aes(x = cohort_label, y = y_pred, fill = dot_col)) +
      geom_jitter(size = 1.2, width = 0.2, shape = point_shape, color = "black", stroke = 0.1) +
      scale_fill_manual(values = dot_col_manual) +
      geom_hline(yintercept = target_99threshold, linetype = "dashed", color = "#161414", linewidth = 0.2) +
      geom_hline(yintercept = target_95threshold, linetype = "dashed", color = "#161414", linewidth = 0.2) +
      theme_classic() +
      theme_NPJ() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Cohort", y = "Score (LR)") +
      # add target_99threshold to y axis breaks and labels
      scale_y_continuous(breaks = c(seq(0, 1, 0.25))) +
      # add a second y axis to the right
      scale_y_continuous(sec.axis = sec_axis(~., breaks = c(0.85, 0.5), labels = c(y_label_for_99spec, y_label_for_95spec))) +
      theme(panel.spacing = unit(0.1, "cm")) +
      theme(strip.text = element_text(size = 6)) +
      theme(axis.text = element_text(size = 6)) +
      theme(axis.title = element_text(size = 6)) +
      theme(legend.position = "none") +
      theme(legend.text = element_text(size = 6)) +
      # lengend title text size = 6
      theme(legend.title = element_text(size = 6)) +
      # set legend title to "Predicted result"
      labs(color = "Predicted result") +
      # remove x anxis title
      theme(axis.title.x = element_blank(), axis.ticks.y.right = element_blank())
  }

  if (target_ichorcna_strat == "(0.1, 1]") {
    p <- ggplot(data = data, aes(x = cohort_label_perfect, y = y_pred, fill = dot_col)) +
      geom_jitter(size = 1.2, width = 0.2, shape = point_shape, color = "black", stroke = 0.1) +
      scale_fill_manual(values = dot_col_manual) +
      geom_hline(yintercept = target_99threshold, linetype = "dashed", color = "#161414", linewidth = 0.2) +
      geom_hline(yintercept = target_95threshold, linetype = "dashed", color = "#161414", linewidth = 0.2) +
      theme_classic() +
      theme_NPJ() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Cohort", y = "Score (LR)") +
      # add target_99threshold to y axis breaks and labels
      scale_y_continuous(breaks = c(seq(0, 1, 0.25))) +
      # add a second y axis to the right
      scale_y_continuous(sec.axis = sec_axis(~., breaks = c(0.5), labels = c("Perfect classifier"))) +
      theme(panel.spacing = unit(0.1, "cm")) +
      theme(strip.text = element_text(size = 6)) +
      theme(axis.text = element_text(size = 6)) +
      theme(axis.title = element_text(size = 6)) +
      theme(legend.position = "none") +
      theme(legend.text = element_text(size = 6)) +
      # lengend title text size = 6
      theme(legend.title = element_text(size = 6)) +
      # set legend title to "Predicted result"
      labs(color = "Predicted result") +
      # remove x anxis title
      theme(axis.title.x = element_blank(), axis.ticks.y.right = element_blank())
  }


  plot_filename <- plot_file
  # get the dir of plot_filename, make the dir if not exists
  plot_dir <- dirname(plot_filename)
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  ggsave(plot_filename,
    p,
    width = 180,
    height = 50,
    unit = "mm",
    create.dir = TRUE
  )

  message(plot_filename)
}

# map the function to the df
pmap_df(df, score_plot, thresholds = thresholds, ans2 = ans2, point_shape = 21)


df_to_save %>%
  write_csv(indtest_csv)

###############################################################################
# SHAP
###############################################################################

# # copy shap from shap_files to /home/nrlab/wang04/ulyses/1_nested_cv_plots/lr_SHAP
# shap_file1 <- file.path(
#   wd,
#   "0-0.03_pairwise",
#   "lr_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0-0.03_pairwise_roc_curve_csv",
#   "SHAP.average_SHAP_summary_plot.dot.pdf"
# )

# shap_file2 <- file.path(
#   wd,
#   "0.03-0.1_pairwise",
#   "lr_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0.03-0.1_pairwise_roc_curve_csv",
#   "SHAP.average_SHAP_summary_plot.dot.pdf"
# )

# shap_file3 <- file.path(
#   wd,
#   "0.1-1_pairwise",
#   "lr_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0.1-1_pairwise_roc_curve_csv",
#   "SHAP.average_SHAP_summary_plot.dot.pdf"
# )
# shap_file4 <- file.path(
#   wd,
#   "all",
#   "lr_input.feat.cnv_sd_length_ctRatio_slRatio.tf.all_roc_curve_csv",
#   "SHAP.average_SHAP_summary_plot.dot.pdf"
# )
# shap_files <- c(shap_file1, shap_file2, shap_file3, shap_file4)

# # # stop if any of the shap_files does not exist
# if (!all(file.exists(shap_files))) {
#   stop("Some of the shap_files do not exist")
# }

# # prepend '0-0.03' to the file name
# shap_file1_des <- gsub("SHAP.average_SHAP_summary_plot.dot.pdf", "0-0.03_SHAP.average_SHAP_summary_plot.dot.pdf", shap_file1)
# # prepend '0.03-0.1' to the file name
# shap_file2_des <- gsub("SHAP.average_SHAP_summary_plot.dot.pdf", "0.03-0.1_SHAP.average_SHAP_summary_plot.dot.pdf", shap_file2)
# # prepend '0.1-1' to the file name
# shap_file3_des <- gsub("SHAP.average_SHAP_summary_plot.dot.pdf", "0.1-1_SHAP.average_SHAP_summary_plot.dot.pdf", shap_file3)
# # prepend 'all' to the file name
# shap_file4_des <- gsub("SHAP.average_SHAP_summary_plot.dot.pdf", "all_SHAP.average_SHAP_summary_plot.dot.pdf", shap_file4)
# # rename the shap_file1 as shap_file1_des
# file.copy(shap_file1, shap_file1_des)
# # rename the shap_file2 as shap_file2_des
# file.copy(shap_file2, shap_file2_des)
# # rename the shap_file3 as shap_file3_des
# file.copy(shap_file3, shap_file3_des)
# # rename the shap_file4 as shap_file4_des
# file.copy(shap_file4, shap_file4_des)


# shap_files_new <- c(shap_file1_des, shap_file2_des, shap_file3_des, shap_file4_des)

# # stop if any of the shap_files_new does not exist
# if (!all(file.exists(shap_files_new))) {
#   stop("Some of the shap_files_new do not exist")
# }

# shap_dir <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/lr_SHAP"
# if (!dir.exists(shap_dir)) {
#   dir.create(shap_dir, recursive = TRUE)
# }
# # copy shap_files to shap_dir
# file.copy(from = shap_files_new, to = shap_dir, overwrite = TRUE)
# message("Finished copying SHAP files to: ", shap_dir)
