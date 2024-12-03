# library(tidyverse)
# library(MultiAssayExperiment)
# library(patchwork)
# library(ggdist)
# library(ggpubr)
# library(gghighlight)
# library(ggmagnify)
# library(plotly)

if (!require("pacman")) install.packages("pacman")

# Load multiple packages
pacman::p_load(tidyverse,
ggupset,
ggplot2,
argparse,
kableExtra,
gridExtra,
grid,
patchwork, 
MultiAssayExperiment,
nord,
plotly,
ggmagnify,
gghighlight,
ggdist,
htmlwidgets,
ggpubr)

source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")
source("/home/nrlab/wang04/ulyses/plot_mae/viz_hyperparmeters.R")

# parameters
timepoints <- 1
plot_file_dir <- "/home/nrlab/wang04/ulyses/plot_mae/plots/ichorcna_tf/"

# mae <- readRDS(mae_file_latest)
# coldata_tibble <- colData(mae) %>% as_tibble()
coldata_tibble <- mae_col_data
raw_coldata_tibble <- read_csv(meta_csv_latest)
# colData_file <- "/home/nrlab/wang04/ulyses/meta_data/colData.csv"
# read in the colData
# coldata_tibble <- read_csv(colData_file)

################################################################################
# stage data clean
################################################################################

dt_stage <- coldata_tibble %>%
  # dplyr::filter(author %in% authors_hyper) %>%
  # dplyr::filter(!is.na(ichorcna_tf)) %>%
  # filter out is.na(stage)
  # dplyr::filter(!is.na(stage)) %>%
  dplyr::mutate(ichorcna_tf = as.numeric(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = round(ichorcna_tf, 3)) %>%
  dplyr::select(author, ichorcna_tf, stage)

# homogenize stage terms
# dt_stage$stage <- dt_stage$stage %>%
#   dplyr::case_match(
#     "0" ~ "Others",
#     "X" ~ "Others",
#     "IA" ~ "I",
#     "IB" ~ "I",
#     "IIA" ~ "II",
#     "IIB" ~ "II",
#     "IIIA" ~ "III",
#     "IIIB" ~ "III",
#     "IIIC" ~ "III",
#     "IVA" ~ "IV",
#     "IVB" ~ "IV",
#     .default = dt_stage$stage
#   )

################################################################################
# ichorcna_tf
################################################################################

max_y <- max(dt_stage$ichorcna_tf)
p_stage <- ggplot(dt_stage, aes(stage, as.numeric(ichorcna_tf))) +
  # ggdist::stat_halfeye(
  #   adjust = .5,
  #   width = .6,
  #   .width = 0,
  #   justification = -.3
  # ) +
  geom_violin(fill = "lightgrey", color = "lightgrey", width = 1, linewidth = 0) +
  geom_point(
    aes(color = author),
    # color ="black",
    size = 0.035,
    alpha = .5,
    # shape =21,
    # stroke = 0.01,
    position = position_jitter(seed = 1, width = .1)
  ) +
  geom_boxplot(
    width = .45,
    outlier.shape = NA,
    width = 0.15,
    linewidth = 0.2,
    alpha = 0
  ) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "lightgrey", linewidth = 0.2) +
  labs(x = "Stage", y = "ichorCNA Tumor Fraction") +
  scale_y_continuous(
    breaks = seq(0.0, 1.0, 0.2),
    # labels = as.character(seq(0.0, 1.0, 0.1)),
    limits = c(0.0, max_y)
  ) +
  theme_classic() +
  # rotate x axis 45
  theme(axis.text = element_text(size = 5.5)) +
  # axis title size to 6
  theme(axis.title = element_text(size = 7)) +
  # strip text size to 6
  theme(strip.text = element_text(size = 5.5)) +
  # strip background rect boarder linewidth to 0.1
  theme(strip.background = element_rect(linewidth = 0, fill = "lightgrey", color = NA)) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  # hide legend
  theme(legend.position = "none") +
  scale_color_manual(values = author_color_values_hyper[unique(dt_stage$author) |> sort()]) +
  # facet by author
  facet_wrap(~author, ncol = 3, scales = "free_y")


# zoom-in to y between 0 and  0.1
y_breaks <- c(0.0, 0.1, 0.03, 0.2, 0.3, 0.4, 0.5)
p_stage_zoom <- p_stage +
  scale_y_continuous(
    breaks = y_breaks,
    labels = as.character(y_breaks),
    limits = c(0.0, 0.1)
  )


plotfile <- file.path(plot_file_dir, "all_stage_tf.pdf")
ggsave(plotfile, plot = p_stage, width = 9, height = 10, units = "cm")
message("Saved plot to ", plotfile)

plotfile <- file.path(plot_file_dir, "all_stage_tf_zoom.pdf")
ggsave(plotfile, plot = p_stage_zoom, width = 9, height = 10, units = "cm")
message("Saved plot to ", plotfile)




################################################################################
# stage bar plot
################################################################################

p_stage_freq <- ggplot(dt_stage, aes(stage, fill = stage)) +
  geom_bar() +
  # add text
  geom_text(
    aes(label = ..count..),
    stat = "count",
    vjust = -0.3,
    size = 2
  ) +
  labs(x = "Stage", y = "Count (n)") +
  scale_y_continuous(expand = c(0.22, 0)) +
  # fill color manual to stage_color_hyper
  scale_fill_manual(values = stage_color_hyper) +
  theme_classic() +
  # make x axis labels smaller
  theme(axis.text = element_text(size = 5.5)) +
  # axis title size to 6
  theme(axis.title = element_text(size = 7)) +
  # strip text size to 6
  theme(strip.text = element_text(size = 5.5)) +
  # strip background rect boarder linewidth to 0.1
  theme(strip.background = element_rect(linewidth = 0, fill = "lightgrey", color = NA)) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(strip.background = element_rect(linewidth = 0, fill = "lightgrey", color = NA)) +
  # hide legend
  theme(legend.position = "none") +
  # strip background color to white
  # theme(strip.background = element_blank()) +
  # scale_fill_manual(values = author_color_values_hyper[unique(dt_stage$author) |> sort()]) +
  facet_wrap(~author, ncol = 3, scales = "free_y")


plotfile <- file.path(plot_file_dir, "stage_bar_facet_by_author.pdf")
ggsave(plotfile, plot = p_stage_freq, width = 9, height = 10, units = "cm")
message("Saved plot to ", plotfile)

################################################################################
# ichorcna_tf vs stage
################################################################################

dt_stage <- dt_stage %>%
  group_by(stage) %>%
  mutate(count = n()) %>%
  mutate(x_label = case_when(
    stage == "I" ~ paste0(stage, "\n", "(n=", count, ")"),
    TRUE ~ paste0(stage, "\n", "(", count, ")")
  ))%>%
  ungroup() %>%
  arrange(stage)

x_levels <- dt_stage$x_label |> unique()

dt_stage$x_label <- factor(dt_stage$x_label, levels = x_levels)

# ggplot y is ichorcna_tf, x is stage
p_stage_tf <- ggplot(dt_stage, aes(x_label, as.numeric(ichorcna_tf))) +
  geom_violin(fill = "lightgrey", color = "lightgrey", width = 1, linewidth = 0) +
  geom_jitter(
    aes(fill = author),
    size = 0.2,
    color = "black",
    alpha = .5,
    stroke = 0.01,
    shape = 21,
    width = 0.2
  ) +
  geom_boxplot(
    width = .45,
    outlier.shape = NA,
    width = 0.15,
    linewidth = 0.2,
    alpha = 0
  ) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "darkgrey", linewidth = 0.2) +
  labs(x = "Stage", y = "ichorCNA Tumor Fraction") +
  scale_y_continuous(
    breaks = c(0.03, 0.1, seq(0.0, 1.0, 0.2)),
    # labels = c("0.03", seq(0.0, 1.0, 0.2)),
    limits = c(0.0, max_y)
  ) +
  theme_classic() +
  theme(axis.text = element_text(size = 5.5)) +
  theme(axis.title = element_text(size = 7)) +
  theme(strip.text = element_text(size = 7)) +
  theme(strip.background = element_rect(linewidth = 0, fill = "lightgrey", color = NA)) +
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(legend.position = "none") +
  scale_fill_manual(values = author_color_values_hyper[unique(dt_stage$author) |> sort()])

y_breaks <- c(0.0, 0.03, 0.1)
p_stage_tf_zoom <- p_stage_tf +
  scale_y_continuous(
    breaks = y_breaks,
    # labels = as.character(y_breaks),
    limits = c(0.0, 0.1)
  )

p_stage_tf_combined <- p_stage_tf +
  (p_stage_tf_zoom +
  #remove y axis title
  theme(axis.title.y = element_blank())
  )

plotfile <- file.path(plot_file_dir, "fig_stage_tf.pdf")
ggsave(plotfile, plot = p_stage_tf, width = 6.5, height = 5, units = "cm")
message("Saved plot to ", plotfile)

plotfile_zoom <- file.path(plot_file_dir, "fig_stage_tf_zoom.pdf")
ggsave(plotfile_zoom, plot = p_stage_tf_zoom, width = 6.5, height = 5, units = "cm")
message("Saved plot to ", plotfile)


plotfile_combined <- file.path(plot_file_dir, "fig_stage_tf_raw_and_zoom_comb.pdf")
ggsave(plotfile_combined, plot = p_stage_tf_combined, width = 11, height = 5, units = "cm")
message("Saved plot to ", plotfile_combined)


################################################################################
# library kit
################################################################################

dt_librarykit <- coldata_tibble %>%
  dplyr::filter(author %in% authors_hyper) %>%
  # dplyr::filter(!is.na(library_kit)) %>%
  dplyr::select(author, library_kit)

# homogenize stage terms



# plot frequency of each stage as bar plot

p_library_kit_freq <- ggplot(dt_librarykit, aes(library_kit)) +
  geom_bar(aes(fill = library_kit)) +
  labs(y = "Frequency (n)") +
  # remove x lab
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  # make x axis labels smaller
  theme(axis.text.x = element_text(size = 6)) +
  # scale_fill_manual(values = author_color_values_hyper[unique(dt_librarykit$author) |> sort()]) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  facet_wrap(~author, ncol = 2, scales = "free_y")


plotfile <- file.path(plot_file_dir, "library_kit_freq.pdf")
ggsave(plotfile, plot = p_library_kit_freq, width = 19, height = 18, units = "cm")
message("Saved plot to ", plotfile)




################################################################################
# extraction kit
################################################################################

dt_extractionkit <- coldata_tibble %>%
  dplyr::filter(author %in% authors_hyper) %>%
  dplyr::select(author, extraction_kit)

# homogenize stage terms



# plot frequency of each stage as bar plot

p_extraction_kit_freq <- ggplot(dt_extractionkit, aes(extraction_kit)) +
  geom_bar(aes(fill = extraction_kit)) +
  labs(y = "Frequency (n)") +
  # remove x lab
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  # make x axis labels smaller
  theme(axis.text.x = element_text(size = 8)) +
  # scale_fill_manual(values = author_color_values_hyper[unique(dt_extractionkit$author) |> sort()]) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  # remove legend
  # theme(legend.position = "none") +
  facet_wrap(~author, ncol = 2, scales = "free_y")


plotfile <- file.path(plot_file_dir, "extraction_kit_freq.pdf")
ggsave(plotfile, plot = p_extraction_kit_freq, width = 18, height = 18, units = "cm")
message("Saved plot to ", plotfile)



################################################################################
# seq_platform
################################################################################

dt_sequencer <- coldata_tibble %>%
  dplyr::filter(author %in% authors_hyper) %>%
  dplyr::select(author, seq_platform)

# homogenize stage terms



# plot frequency of each stage as bar plot

p_sequencer_freq <- ggplot(dt_sequencer, aes(seq_platform)) +
  geom_bar(aes(fill = seq_platform)) +
  labs(y = "Frequency (n)") +
  # remove x lab
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  # make x axis labels smaller
  theme(axis.text.x = element_text(size = 8)) +
  # scale_fill_manual(values = author_color_values_hyper[unique(dt_extractionkit$author) |> sort()]) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  # remove legend
  # theme(legend.position = "none") +
  facet_wrap(~author, ncol = 2, scales = "free_y")


plotfile <- file.path(plot_file_dir, "seq_platform_freq.pdf")
ggsave(plotfile, plot = p_sequencer_freq, width = 20, height = 18, units = "cm")
message("Saved plot to ", plotfile)


################################################################################
# gender
################################################################################

dt_gender <- coldata_tibble %>%
  dplyr::filter(author %in% authors_hyper) %>%
  dplyr::select(author, gender)

# homogenize stage terms

dt_gender$gender <- dt_gender$gender %>%
  dplyr::case_match(
    "N/A" ~ NA,
    "F" ~ "female",
    "M" ~ "male",
    .default = dt_gender$gender
  )


# plot frequency of each stage as bar plot

p_gender_freq <- ggplot(dt_gender, aes(gender)) +
  geom_bar(aes(fill = gender)) +
  labs(y = "Frequency (n)") +
  # remove x lab
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  # make x axis labels smaller
  theme(axis.text.x = element_text(size = 8)) +
  # scale_fill_manual(values = author_color_values_hyper[unique(dt_extractionkit$author) |> sort()]) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  # remove legend
  # theme(legend.position = "none") +
  facet_wrap(~author, ncol = 2, scales = "free_y")


plotfile <- file.path(plot_file_dir, "gender_freq.pdf")
ggsave(plotfile, plot = p_gender_freq, width = 20, height = 18, units = "cm")
message("Saved plot to ", plotfile)

################################################################################
# age
################################################################################

dt_age <- coldata_tibble %>%
  dplyr::filter(author %in% authors_hyper) %>%
  dplyr::select(author, age)



# plot frequency of each stage as bar plot

p_age_freq <- ggplot(dt_age, aes(age)) +
  geom_bar() +
  labs(y = "Frequency (n)") +
  # remove x lab
  theme_classic() +
  theme(axis.title.x = element_blank()) +
  # make x axis labels smaller
  theme(axis.text.x = element_text(size = 8)) +
  # scale_fill_manual(values = author_color_values_hyper[unique(dt_extractionkit$author) |> sort()]) +
  # rotate x axis labels 45 degree
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  # remove legend
  # theme(legend.position = "none") +
  facet_wrap(~author, ncol = 2, scales = "free_y")


plotfile <- file.path(plot_file_dir, "age_freq.pdf")
ggsave(plotfile, plot = p_age_freq, width = 20, height = 18, units = "cm")
message("Saved plot to ", plotfile)

################################################################################
# box plot of ichorcna_tf stratified by author
################################################################################

authors <- authors_hyper
# cancer
dt_cancer <- coldata_tibble %>%
  dplyr::filter(author %in% authors) %>%
  # filter out non-cancer disease cohort
  dplyr::filter(!author %in% noncancer_disease_hyper) %>%
  dplyr::filter(!is.na(ichorcna_tf)) %>%
  # filter out "Healthy" cohort
  dplyr::filter(cohort != "Healthy") %>%
  dplyr::mutate(ichorcna_tf = as.numeric(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = round(ichorcna_tf, 3)) %>%
  dplyr::select(author, ichorcna_tf, stage)



p <- ggplot(dt_cancer, aes(forcats::fct_reorder(author, ichorcna_tf), ichorcna_tf)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .25,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = author),
    size = 0.6,
    alpha = .4,
    position = position_jitter(seed = 1, width = .1)
  ) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "red") +
  labs(x = "Dataset", y = "ichorCNA Tumor Fraction") +
  # set color manually
  scale_color_manual(values = author_color_values_hyper[unique(dt$author) |> sort()]) +
  scale_y_continuous(
    breaks = seq(0.0, 1.0, 0.1),
    labels = as.character(seq(0.0, 1.0, 0.1)),
    limits = c(0.0, 1.0)
  ) +
  # hide legend
  guides(color = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# zoom-in to y between 0 and  0.5
y_breaks <- c(0.1, 0.03, 0.2, 0.3, 0.4, 0.5)
p_zoom <- p +
  scale_y_continuous(
    breaks = y_breaks,
    labels = as.character(y_breaks),
    limits = c(0.0, 0.5)
  )
plotfile <- file.path(plot_file_dir, "author_cancer_tf.pdf")
ggsave(plotfile, plot = p, width = 18, height = 8, units = "cm")
message("Saved plot to ", plotfile)

plotfile <- file.path(plot_file_dir, "author_cancer_tf_zoom.pdf")
ggsave(plotfile, plot = p_zoom, width = 18, height = 8, units = "cm")
message("Saved plot to ", plotfile)


# healthy
dt_healthy <- coldata_tibble %>%
  dplyr::filter(author %in% authors) %>%
  # filter out non-cancer disease cohort
  dplyr::filter(!is.na(ichorcna_tf)) %>%
  # filter out "Healthy" cohort
  dplyr::filter(cohort == "Healthy") %>%
  dplyr::mutate(ichorcna_tf = as.numeric(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = round(ichorcna_tf, 3)) %>%
  dplyr::select(author, ichorcna_tf, stage)



p <- ggplot(dt_healthy, aes(forcats::fct_reorder(author, ichorcna_tf), ichorcna_tf)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .25,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = author),
    size = 0.6,
    alpha = .4,
    position = position_jitter(seed = 1, width = .1)
  ) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "red") +
  labs(x = "Dataset", y = "ichorCNA Tumor Fraction") +
  # set color manually
  scale_color_manual(values = author_color_values_hyper[unique(dt$author) |> sort()]) +
  scale_y_continuous(
    breaks = seq(0.0, 1.0, 0.1),
    labels = as.character(seq(0.0, 1.0, 0.1)),
    limits = c(0.0, 1.0)
  ) +
  # hide legend
  guides(color = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# zoom-in to y between 0 and  0.5
y_breaks <- c(0, 0.01, 0.02, 0.03, 0.04)
p_zoom <- p +
  scale_y_continuous(
    breaks = y_breaks,
    labels = as.character(y_breaks),
    limits = c(0.0, 0.04)
  )
plotfile <- file.path(plot_file_dir, "author_tf_healthy.pdf")
ggsave(plotfile, plot = p, width = 18, height = 8, units = "cm")
message("Saved plot to ", plotfile)

plotfile <- file.path(plot_file_dir, "author_tf_healthy_zoom.pdf")
ggsave(plotfile, plot = p_zoom, width = 18, height = 8, units = "cm")
message("Saved plot to ", plotfile)


# NCD

dt_NCD <- coldata_tibble %>%
  dplyr::filter(author %in% authors) %>%
  # filter out non-cancer disease cohort
  dplyr::filter(!is.na(ichorcna_tf)) %>%
  # filter out "Healthy" cohort
  dplyr::filter(cohort %in% noncancer_disease_hyper) %>%
  dplyr::mutate(ichorcna_tf = as.numeric(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = round(ichorcna_tf, 3)) %>%
  dplyr::select(author, ichorcna_tf, stage)



p <- ggplot(dt_NCD, aes(forcats::fct_reorder(author, ichorcna_tf), ichorcna_tf)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA
  ) +
  geom_boxplot(
    width = .25,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = author),
    size = 0.6,
    alpha = .4,
    position = position_jitter(seed = 1, width = .1)
  ) +
  geom_hline(yintercept = 0.03, linetype = "dashed", color = "red") +
  labs(x = "Dataset", y = "ichorCNA Tumor Fraction") +
  # set color manually
  scale_color_manual(values = author_color_values_hyper[unique(dt$author) |> sort()]) +
  scale_y_continuous(
    breaks = seq(0.0, 1.0, 0.1),
    labels = as.character(seq(0.0, 1.0, 0.1)),
    limits = c(0.0, 1.0)
  ) +
  # hide legend
  guides(color = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# zoom-in to y between 0 and  0.5
y_breaks <- c(0, 0.01, 0.02, 0.03, 0.04)
p_zoom <- p +
  scale_y_continuous(
    breaks = y_breaks,
    labels = as.character(y_breaks),
    limits = c(0.0, 0.04)
  )
plotfile <- file.path(plot_file_dir, "author_tf_NCD.pdf")
ggsave(plotfile, plot = p, width = 18, height = 8, units = "cm")
message("Saved plot to ", plotfile)

plotfile <- file.path(plot_file_dir, "author_tf_NCD_zoom.pdf")
ggsave(plotfile, plot = p_zoom, width = 18, height = 8, units = "cm")
message("Saved plot to ", plotfile)

################################################################################
# compare all clinical_tf with ichorcna_tf
################################################################################

x_breaks <- c(0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
y_breaks <- seq(0, 1.0, 0.1)

coldata_tibble <- as_tibble(colData(mae), rownames = "bam_id")

dt_tf <- coldata_tibble %>%
  dplyr::filter(author %in% authors_hyper) %>%
  dplyr::filter(!is.na(ichorcna_tf)) %>%
  # only keep clinical_tf not NA
  dplyr::filter(!is.na(clinical_tf)) %>%
  dplyr::mutate(ichorcna_tf = as.numeric(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = round(ichorcna_tf, 3)) %>%
  dplyr::mutate(clinical_tf = round(clinical_tf, 3)) %>%
  dplyr::select(author, ichorcna_tf, clinical_tf, bam_id)


p_tf <- ggplot(dt_tf, aes(x = ichorcna_tf, y = clinical_tf)) +
  geom_point(
    color = "black",
    fill = "#e2dede",
    shape = 21
  ) +
  geom_point(
    data = dt_tf %>% filter(ichorcna_tf < 0.03),
    color = "black",
    fill = "salmon",
    shape = 21
  ) +
  geom_point(
    data = dt_tf %>% filter(clinical_tf == 0),
    color = "black",
    fill = "orange",
    alpha = 0.5,
    shape = 21
  ) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "black") +
  labs(x = "ichorCNA Tumor Fraction", y = "Clinical Tumor Fraction") +
  scale_x_continuous(
    breaks = x_breaks,
    labels = as.character(x_breaks),
    limits = c(0.0, 1.0)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = as.character(y_breaks),
    limits = c(0.0, 1.0)
  ) +
  theme_classic() +
  # rotate x axis labels
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggpubr::stat_cor(method = "pearson", label.x = 0.7, label.y = 0.9, cor.coef.name = "R") +
  # facet by author
  facet_wrap(~author, ncol = 1, scales = "free_y")

outfile <- file.path(plot_file_dir, "all_tf_comparison.pdf")
ggsave(outfile, plot = p_tf, width = 18, height = 18, units = "cm")
message("Saved plot to ", outfile)



################################################################################
# compare pp et al  clinical_tf with ichorcna_tf
################################################################################
authors <- c("P.P. et al")

# set y breaks
x_breaks <- c(0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
y_breaks <- seq(0, 1.0, 0.1)

coldata_tibble <- as_tibble(colData(mae), rownames = "bam_id")

dt <- coldata_tibble %>%
  dplyr::filter(author %in% authors) %>%
  dplyr::filter(!is.na(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = as.numeric(ichorcna_tf)) %>%
  dplyr::mutate(ichorcna_tf = round(ichorcna_tf, 3)) %>%
  dplyr::mutate(clinical_tf = round(clinical_tf, 3)) %>%
  dplyr::select(ichorcna_tf, clinical_tf, bam_id)

p <- ggplot(dt, aes(x = ichorcna_tf, y = clinical_tf)) +
  geom_point(
    color = "black",
    fill = "#e2dede",
    shape = 21
  ) +
  geom_point(
    data = dt %>% filter(ichorcna_tf < 0.03),
    color = "black",
    fill = "salmon",
    shape = 21
  ) +
  geom_point(
    data = dt %>% filter(clinical_tf == 0),
    color = "black",
    fill = "orange",
    alpha = 0.5,
    shape = 21
  ) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "black") +
  labs(x = "ichorCNA Tumor Fraction", y = "Clinical Tumor Fraction") +
  scale_x_continuous(
    breaks = x_breaks,
    labels = as.character(x_breaks),
    limits = c(0.0, 1.0)
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    labels = as.character(y_breaks),
    limits = c(0.0, 1.0)
  ) +
  theme_classic() +
  # rotate x axis labels
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggpubr::stat_cor(method = "pearson", label.x = 0.7, label.y = 0.9, cor.coef.name = "R") +
  geom_magnify(
    from = c(0, 0, 0.1, 0.1),
    to = c(0.41, 0, 1, 0.4),
    shadow = TRUE,
    proj.linetype = 3,
    colour = "lightgray",
    axes = "xy"
  )

outfile <- file.path(plot_file_dir, "pp_et_al_tf_comparison.pdf")
ggsave(outfile, plot = p, width = 8, height = 6, units = "in")
message("Saved plot to ", outfile)



################################################################################
# density plot of n_frag stratified by data source
################################################################################
library(tidyverse)
library(ggplot2)
library(gghighlight)
timepoints <- c(1)

# colData_file <- "/home/nrlab/wang04/ulyses/meta_data/colData.csv"
# coldata_tibble <- read_csv(colData_file)

author <- unique(coldata_tibble$author)

target_author <- authors_hyper
target_depth <- seq(0.1, 5, 0.1)
# generate all combinations of author and target_depth
author_target_depth <- crossing(target_author, target_depth)

coldata_subset <- coldata_tibble %>%
  dplyr::filter(.data$sample_type != "urine") %>%
  dplyr::filter(.data$author %in% target_author)
# dplyr::filter(is.na(timepoint) | timepoint %in% timepoints)


filter_depth <- function(target_author, target_depth) {
  result <- coldata_subset %>%
    dplyr::filter(.data$author == !!target_author) %>%
    dplyr::filter(.data$seq_depth >= as.numeric(!!target_depth)) %>%
    nrow()
  return(result)
}

rslt <- purrr::pmap(author_target_depth, filter_depth) %>%
  unlist()

author_target_depth <- author_target_depth %>%
  mutate(n = rslt)

p1 <- ggplot(coldata_subset, aes(x = seq_depth)) +
  geom_histogram(fill = "#796e6e", binwidth = 0.1) +
  labs(x = "Sequencing Depth", y = "N Samples") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  theme_classic() +
  coord_cartesian(xlim = c(0, 50)) +
  facet_wrap(~author, ncol = 3, scales = "free_y")

p2 <- ggplot(coldata_subset, aes(x = seq_depth)) +
  geom_density(fill = "#e3c5c5", alpha = 0.5) +
  labs(x = "Sequencing Depth", y = "Density") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  theme_classic() +
  coord_cartesian(xlim = c(0, 50)) +
  facet_wrap(~author, ncol = 3, scales = "free_y")

p3 <- ggplot(author_target_depth, aes(x = as.factor(target_depth), y = n)) +
  geom_bar(stat = "identity", fill = "#b9a0a0") +
  # set the bar where x = 0.1 to orange
  geom_bar(
    data = author_target_depth %>% filter(target_depth == 0.1),
    stat = "identity", fill = "#161517"
  ) +
  labs(x = "Coverage Cutoff (X)", y = "N Samples") +
  theme_classic() +
  # only show the x axis text for 1, 2 , 3, 4, 5
  scale_x_discrete(breaks = c(1, 2, 3, 4, 5)) +
  # rotate the x axis labels -90 degree
  theme(axis.text.x = element_text(angle = 10, hjust = 1)) +
  facet_wrap(~target_author, ncol = 3, scales = "free_y")

# stacked bar chart, fill by author
p4 <- ggplot(author_target_depth, aes(x = as.factor(target_depth), y = n, fill = target_author)) +
  geom_bar(stat = "identity") +
  labs(x = "Coverage Cutoff (X)", y = "N Samples") +
  theme_classic() +
  # set color manually
  scale_fill_manual(values = author_color_values_hyper[unique(author_target_depth$target_author) |> sort()]) +
  # only show the x axis text for 1, 2 , 3, 4, 5
  scale_x_discrete(breaks = c(1, 2, 3, 4, 5)) +
  # rotate the x axis labels -90 degree
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

p1_file <- file.path(plot_file_dir, "seq_depth_histogram.pdf")
ggsave(p1_file, plot = p1, width = 18, height = 8, units = "cm")
message("Saved plot to ", p1_file)

p2_file <- file.path(plot_file_dir, "seq_depth_density.pdf")
ggsave(p2_file, plot = p2, width = 18, height = 8, units = "cm")
message("Saved plot to ", p2_file)

p3_file <- file.path(plot_file_dir, "seq_depth_bar.pdf")
ggsave(p3_file, plot = p3, width = 18, height = 8, units = "cm")
message("Saved plot to ", p3_file)

p4_file <- file.path(plot_file_dir, "seq_depth_stacked_bar.pdf")
ggsave(p4_file, plot = p4, width = 18, height = 8, units = "cm")
message("Saved plot to ", p4_file)


###############################################################################
# plot density of n timepoints stratified by author
###############################################################################
library(tidyverse)
library(ggplot2)
library(gghighlight)
timepoints <- c(1)


# colData_file <- "/home/nrlab/wang04/ulyses/meta_data/colData.csv"
# coldata_tibble <- read_csv(colData_file)

authors <- authors_hyper
author <- unique(coldata_tibble$author)

target_author <- authors_hyper

coldata_subset <- raw_coldata_tibble %>%
  dplyr::filter(.data$sample_type %in% sample_type_hyper) %>%
  dplyr::filter(.data$author %in% target_author) %>%
  # dplyr::filter(.data$seq_depth >= 0.1) %>%
  dplyr::select(author, timepoint, patient_id)

ans <- coldata_subset %>%
  group_by(author, patient_id) %>%
  summarise(n_timepoints = n()) %>%
  ungroup()

# set n_timepoints to character

p <- ggplot(ans, aes(x = n_timepoints)) +
  geom_histogram(aes(fill = author), binwidth = 1) +
  labs(x = "N Timepoints", y = "N Samples") +
  theme_classic() +
  scale_color_manual(values = author_color_values_hyper[unique(ans$author) |> sort()]) +
  facet_wrap(~author, ncol = 3, scales = "free_y")



# stacked bar chart, fill by author
p_stacked <- ggplot(ans, aes(x = n_timepoints, fill = author)) +
  geom_bar(aes(fill = author)) +
  labs(x = "N Timepoints", y = "N Samples") +
  coord_cartesian(ylim = c(0, 1500)) +
  scale_color_manual(values = author_color_values_hyper[unique(ans$author) |> sort()]) +
  theme_classic()

p_file <- file.path(plot_file_dir, "n_timepoints_histogram.pdf")
ggsave(p_file, plot = p, width = 18, height = 12, units = "cm")
message("Saved plot to ", p_file)

p_stacked_file <- file.path(plot_file_dir, "n_timepoints_stacked_bar.pdf")
ggsave(p_stacked_file, plot = p_stacked, width = 18, height = 12, units = "cm")
message("Saved plot to ", p_stacked_file)
