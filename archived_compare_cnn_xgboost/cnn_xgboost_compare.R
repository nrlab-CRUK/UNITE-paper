
library(ggupset)
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(argparse)
library(kableExtra)
library(gridExtra)
library(grid)
library(nord)


# make parameters using argparse
parser <- argparse::ArgumentParser()
parser$add_argument("--xgboost_wd", type = "character", help = "xgboost working directory", default = "/scratchc/nrlab/wang04/ulyses/xgboost_model_v1")
parser$add_argument("--cnn_wd", type = "character", help = "cnn working directory", default = "/scratchc/nrlab/wang04/ulyses/cnn_model/cnn_strat/finalize/tp1")
parser$add_argument("--which_xgboost_feat", type = "character", help = "which feature to plot", default = "cnv_sd_length_motif")
parser$add_argument("--cnn_file_pattern", type = "character", help = "cnn file pattern", default = "performance.*\\.debug.csv$")
parser$add_argument("--cnn_metrics_to_plot", type = "character", help = "cnn metrics to plot", default = c("accuracy", "AUROC"))
parser$add_argument("--ichorcna_strat_to_plot", type = "character", help = "ichorcna strata to plot", default = c("(0, 0.01]", "(0.01, 0.03]", "(0.03, 0.1]", "(0.1, 0.2]", "(0.2, 1]", "all"))
parser$add_argument("--parallel_coord_plot_file", type = "character", help = "parallel coordinate plot file", default = file.path("/home/nrlab/wang04/ulyses/compare_cnn_xgboost", "cnn_xgboost_parallel_coord_plot.pdf"))

# add ichorcna_strat_levels, default = c("(0, 0.01]", "(0.01, 0.03]", "(0.03, 0.1]", "(0.1, 0.2]", "(0.2, 1]", "all")
parser$add_argument("--ichorcna_strat_levels", type = "character", help = "ichorcna strata levels", default = c("(0, 0.01]", "(0.01, 0.03]", "(0.03, 0.1]", "(0.1, 0.2]", "(0.2, 1]", "all"))
# add all_hight_color, default = "#c1dde7"
parser$add_argument("--all_hight_color", type = "character", help = "all hight color", default = "#bcdae5")

# add labs_x, default = "ichorCNA TF Strata"
parser$add_argument("--labs_x", type = "character", help = "labs x", default = NULL)
# add labs_y, default = "Mean Performance"
parser$add_argument("--labs_y", type = "character", help = "labs y", default = "Mean Performance")

# add xgboost_cnn_manual_color, default = c("#af8181", "#cd1690")
parser$add_argument("--xgboost_cnn_manual_color", type = "character", help = "xgboost cnn manual color", default = c("#af8181", "#1c7de4"))

# add plot_width, default = 180
parser$add_argument("--plot_width", type = "integer", help = "plot width", default = 100)
# add plot_height, default = 80
parser$add_argument("--plot_height", type = "integer", help = "plot height", default = 50)
# add plot_units, default = "mm"
parser$add_argument("--plot_units", type = "character", help = "plot units", default = "mm")
# add plot_dpi, default = 300
parser$add_argument("--plot_dpi", type = "integer", help = "plot dpi", default = 300)

# get the parameters
args <- parser$parse_args()
xgboost_wd <- args$xgboost_wd
cnn_wd <- args$cnn_wd
which_xgboost_feat <- args$which_xgboost_feat
cnn_file_pattern <- args$cnn_file_pattern
cnn_metrics_to_plot <- args$cnn_metrics_to_plot
ichorcna_strat_to_plot <- args$ichorcna_strat_to_plot
parallel_coord_plot_file <- args$parallel_coord_plot_file
ichorcna_strat_levels <- args$ichorcna_strat_levels
all_hight_color <- args$all_hight_color
labs_x <- args$labs_x
labs_y <- args$labs_y
xgboost_cnn_manual_color <- args$xgboost_cnn_manual_color
plot_width <- args$plot_width
plot_height <- args$plot_height
plot_units <- args$plot_units
plot_dpi <- args$plot_dpi


# get the xgboost table

fl <- list.files(path = xgboost_wd,
                 "repeat\\d+_.*\\.csv", 
		 recursive = TRUE,
                 full.names=TRUE)
fl_read <- lapply(fl, read_csv, show_col_types = FALSE)
names(fl_read) <- fl
fl_read_bind <- bind_rows(fl_read, .id = "id") %>%
  select(-c(.estimator, .config)) %>%
  mutate(ichorcna_strat = case_when(
	stringr::str_detect(id, "0-0.01") ~ '(0, 0.01]',
	stringr::str_detect(id, "0.01-0.03") ~ '(0.01, 0.03]',
	stringr::str_detect(id, "0.03-0.1") ~ '(0.03, 0.1]',
	stringr::str_detect(id, "0.1-0.2") ~ '(0.1, 0.2]',
	stringr::str_detect(id, "0.2-1") ~ '(0.2, 1]',
	stringr::str_detect(id, "all") ~ 'all',
	TRUE ~ NA_character_
  ))
# set order of ichorcna_strat
fl_read_bind$ichorcna_strat <- factor(fl_read_bind$ichorcna_strat, levels = ichorcna_strat_levels)


# plot the feature performance of 'cnv_sd_length_motif'
ans <- fl_read_bind %>%
filter(feat == !!which_xgboost_feat) %>%
select(-c(id, feat))


# summarize as table

xgboost_tibble <- ans %>% 
  group_by(ichorcna_strat, .metric) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  mutate(`to_report` = paste(round(mean, 3), round(sd, 3), sep = "±")) %>%
  select(-c(mean, sd)) %>% 
  #pivot_wider(names_from = `ichorcna_strat`, values_from = to_report) %>%
  rename(metrics = .metric) %>%
  # change 'roc_auc' to 'AUROC'
  mutate(metrics = ifelse(metrics == 'roc_auc', 'AUROC', metrics)) %>%
  mutate(model = "XGBoost")

###############################################################################
# get the cnn table
###############################################################################


fl <- list.files(path = cnn_wd, pattern = cnn_file_pattern, full.names = TRUE)

fl_read <- map(fl, read_csv, show_col_types = FALSE) %>% 
  map(pivot_longer, cols = everything(), names_to = "metrics", values_to = "value") %>%
  setNames(fl) %>%
  bind_rows(.id = "file") %>%
  mutate(`ichorCNA Tumor Fraction` = case_when(
    str_detect(file, "0_0.01") ~ "(0, 0.01]",
    str_detect(file, "0.01_0.03") ~ "(0.01, 0.03]",
    str_detect(file, "0.03_0.1") ~ "(0.03, 0.1]",
    str_detect(file, "0.1_0.2") ~ "(0.1, 0.2]",
    str_detect(file, "0.2_1") ~ "(0.2, 1]",
    str_detect(file, "0_1") ~ "all"
  )) 



fl_read_train_test <- fl_read %>% 
  mutate(train_test = case_when(
    str_detect(metrics, "train") ~ "train",
    str_detect(metrics, "test") ~ "test",
    TRUE ~ "NA"
  )) %>%
  mutate(metrics = str_remove(metrics, "^train_|^test_"))

fl_read_train_test$train_test <- factor(fl_read_train_test$train_test, levels = c("train", "test"))

# filter tibble using metric_to_plot
fl_read_train_test_use <- fl_read_train_test %>% filter(metrics %in% cnn_metrics_to_plot)

cnn_tibble <- fl_read_train_test_use %>% 
  group_by(`ichorCNA Tumor Fraction`, train_test, metrics) %>%
  summarize(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(`to_report` = paste(round(mean, 3), round(sd, 3), sep = "±")) %>%
  select(-c(mean, sd)) %>% 
  rename(ichorcna_strat = `ichorCNA Tumor Fraction`) %>%
  #pivot_wider(names_from = ichorcna_strat , values_from = to_report) %>%
  filter(train_test == "test") %>%
  select(-c(train_test)) %>%
  mutate(model = "CNN")


###############################################################################
# combine the two tables
###############################################################################
combined_tibble <- bind_rows(xgboost_tibble, cnn_tibble) %>%
  mutate(model = factor(model, levels = c("XGBoost", "CNN"))) %>%
  # change 'accuracy' to 'Accuracy'
  mutate(metrics = ifelse(metrics == 'accuracy', 'Accuracy', metrics))

# remove the part after '±' in to_report and set it as numeric
combined_tibble$to_report <- as.numeric(str_remove(combined_tibble$to_report, "±.*"))

# set order of ichorcna_strat
combined_tibble$ichorcna_strat <- factor(combined_tibble$ichorcna_strat, 
				levels = ichorcna_strat_levels)

###############################################################################
# plot as parallel coordinates plot
###############################################################################

plot_data <- combined_tibble %>%
	filter(ichorcna_strat %in% ichorcna_strat_to_plot)

# plot_data without 'all' group
plot_data_no_all <- plot_data %>% filter(ichorcna_strat != "all")

p <- ggplot(plot_data, aes(x = ichorcna_strat, y = to_report, color = model, group = model)) +
  geom_line(data = plot_data_no_all) +
  # Add background color to 'all' group
  geom_rect(
    data = subset(plot_data, ichorcna_strat == "all"),
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
  geom_point() +
  facet_wrap(~metrics) +
  # set y axis as three decimal places
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  #facet_grid(metrics ~ .) +
  scale_color_manual(values = xgboost_cnn_manual_color) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # add legend
  theme(strip.background = element_rect(fill = "#e9e6e6", color = "#e9e6e6")) +
  labs(x = labs_x, y = labs_y) +
  #theme(panel.spacing = unit(0.1, "cm")) +
  theme(strip.text = element_text(size = 5)) +
  theme(axis.text = element_text(size = 5)) +
  # remove legend title
  theme(axis.title = element_text(size = 7)) +
  theme(legend.position = c(.80, .25)) +
  theme(legend.text = element_text(size = 5), legend.title = element_blank()) +
  #theme(legend.key.size = unit(0.5, "cm")) +
  #theme(legend.key.width = unit(0.5, "cm")) +
  theme(legend.key = element_rect(fill = "transparent", color = "transparent")) +
  theme(legend.background = element_rect(fill = "transparent", color = "transparent")) +
  theme(legend.box.background = element_rect(fill = "transparent", color = "transparent")) +
  # legend nearer the plot
  scale_x_discrete(expand = c(0.05, 0))


# save plot
ggsave(
  filename = parallel_coord_plot_file,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = plot_units,
  dpi = plot_dpi)

message("Saved to ", parallel_coord_plot_file)
