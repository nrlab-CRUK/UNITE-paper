library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(ComplexHeatmap)

cnv_plot_saving_dir_base <- "/home/nrlab/wang04/ulyses/plot_mae/plots/cnv"

source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")
source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

# make cnv_plot_saving_dir with iteration_num_hyper
cnv_plot_saving_dir <- file.path(cnv_plot_saving_dir_base, paste0("iteration_", iteration_num_hyper))
if (!dir.exists(cnv_plot_saving_dir)) {
  dir.create(cnv_plot_saving_dir)
}



################################################################################
# cnv
################################################################################
all_cnv <- all_long_filter %>%
  dplyr::filter(assay %in% c("cnv")) %>%
  # remove rows 6p_6
  dplyr::filter(rowname != "6p_6") %>%
  dplyr::mutate(chr_arm = str_extract(rowname, "\\d+[pq]")) %>%
  dplyr::mutate(chr = str_extract(chr_arm, "\\d+")) %>%
  dplyr::mutate(arm = str_extract(chr_arm, "[pq]")) %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)


all_cnv_grouped <- all_cnv %>%
  group_by(ichorcna_tf_strat)

group_keys(all_cnv_grouped)

tf_list <- all_cnv %>%
  group_by(ichorcna_tf_strat) %>%
  group_split()

names(tf_list) <- group_keys(all_cnv_grouped) |> pull(ichorcna_tf_strat)

make_matrix <- function(df, rownames = NULL) {
  my_matrix <- as.matrix(df)
  if (!is.null(rownames)) {
    rownames(my_matrix) <- rownames
  }
  my_matrix
}
f_wider <- function(tb) {
  ans <- tb %>%
    select(primary, x, y) %>%
    # convert to wider format using pivot_wider
    pivot_wider(names_from = primary, values_from = y)

  ans2 <- make_matrix(ans %>% select(-x), pull(ans, x))

  # ans_scaled <- apply(ans2, 2, scale)
  ans2[] <- apply(ans2, 2, scale)
  return(ans2)
}

tf_list_wide <- lapply(tf_list, f_wider)


saveRDS(tf_list_wide, file.path(cnv_plot_saving_dir, "tf_list_wide.rds"))

# Heatmap

strat_plot <- file.path(cnv_plot_saving_dir, "cnv_heatmap_split")
one_off_plot <- file.path(cnv_plot_saving_dir, "cnv_heatmap_all")

# plot as heatmap

source("/home/nrlab/wang04/ulyses/plot_mae/heatmap.R")


h1 <- plot_heatmap(tf_list_wide, 1,
  annotation_name_side = "left",
  show_detailed_sample_annotation = FALSE
)

h2 <- plot_heatmap(tf_list_wide, 2,
  show_heatmap_legend = FALSE,
  show_chr_name = FALSE,
  show_detailed_sample_annotation = FALSE,
  show_annotation_name = FALSE,
  show_point_anno_axis = FALSE
)

h3 <- plot_heatmap(tf_list_wide, 3,
  show_heatmap_legend = FALSE,
  show_chr_name = FALSE,
  show_detailed_sample_annotation = FALSE,
  show_annotation_name = FALSE,
  show_point_anno_axis = FALSE
)

h4 <- plot_heatmap(tf_list_wide, 4,
  show_heatmap_legend = FALSE,
  show_chr_name = FALSE,
  show_detailed_sample_annotation = FALSE,
  show_annotation_name = FALSE,
  show_point_anno_axis = FALSE
)

# h5 <- plot_heatmap(tf_list_wide, 5,
#   show_heatmap_legend = FALSE,
#   show_chr_name = FALSE,
#   show_detailed_sample_annotation = FALSE,
#   show_annotation_name = FALSE,
#   show_point_anno_axis = FALSE
# )


h_combine <- h1 + h2 + h3 + h4
# plot as heatmap

pdf(paste0(strat_plot, ".pdf"), width = 4.88, height = 3.78)
draw(h_combine, ht_gap = unit(0.2, "cm"), merge_legend = TRUE)
dev.off()


# plot all samples at a time ----------------------------------------------


# combine everything
tf_list_wide_all <- purrr::reduce(tf_list_wide, cbind) %>% list()
hall <- plot_heatmap(tf_list_wide_all, 1, annotation_name_side = "left")

# pdf(paste0(one_off_plot, ".pdf"), width = 4.88, height = 3.78)
pdf(paste0(one_off_plot, ".pdf"), width = 18, height = 9.6)
draw(hall)
dev.off()


png(paste0(one_off_plot, ".png"), width = 18, height = 9.6, units = "in", res = 300)
draw(hall)
dev.off()


################################################################################
# prepare cnv input data (segmented)
################################################################################


#############################################################
# segment cnv
#############################################################

all_cnv_seg <- all_long_filter %>%
  dplyr::filter(assay == "cnv") %>%
  dplyr::filter(rowname != "6p_6") %>%
  dplyr::mutate(chr = str_extract(rowname, "\\d+")) %>%
  group_by(primary, chr) %>%
  dplyr::mutate(pos = row_number()) %>%
  ungroup() %>%
  dplyr::group_by(primary)

all_cnv_seg$rowname <- factor(all_cnv_seg$rowname, levels = all_cnv_seg$rowname %>% unique())

all_cnv_seg_list <- group_split(all_cnv_seg)
names(all_cnv_seg_list) <- group_keys(all_cnv_seg) %>% pull()

# seg function from /home/nrlab/wang04/ulyses/plot_mae/heatmap.R
# and already sourced above

all_cnv_seg <- lapply(all_cnv_seg_list, seg) %>%
  bind_rows() %>%
  select(-c(chr, pos)) %>%
  dplyr::filter(assay %in% c("cnv")) %>%
  dplyr::filter(rowname != "6p_6") %>%
  dplyr::mutate(chr_arm = str_extract(rowname, "\\d+[pq]")) %>%
  dplyr::mutate(chr = str_extract(chr_arm, "\\d+")) %>%
  dplyr::mutate(arm = str_extract(chr_arm, "[pq]")) %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)


# IMPORTANT: the levels of bicohort affects how line 351 is set.
all_cnv_seg$bicohort <- factor(all_cnv_seg$bicohort, levels = c("Healthy", "Cancer"))


###

all_cnv_seg_grouped <- all_cnv_seg %>%
  group_by(ichorcna_tf_strat)

group_keys(all_cnv_seg_grouped)

tf_list_seg <- all_cnv_seg %>%
  group_by(ichorcna_tf_strat) %>%
  group_split()

names(tf_list_seg) <- group_keys(all_cnv_seg_grouped) |> pull(ichorcna_tf_strat)

make_matrix <- function(df, rownames = NULL) {
  my_matrix <- as.matrix(df)
  if (!is.null(rownames)) {
    rownames(my_matrix) <- rownames
  }
  my_matrix
}

f_wider_no_scale <- function(tb) {
  ans <- tb %>%
    select(primary, x, y) %>%
    # convert to wider format using pivot_wider
    pivot_wider(names_from = primary, values_from = y)

  ans2 <- make_matrix(ans %>% select(-x), pull(ans, x))

  # ans_scaled <- apply(ans2, 2, scale)
  # ans2[] <- apply(ans2, 2, scale)
  return(ans2)
}

tf_list_seg_wide <- lapply(tf_list_seg, f_wider_no_scale)


saveRDS(tf_list_seg_wide, file.path(cnv_plot_saving_dir, "tf_list_seg_wide.rds"))


# do seg plot
# Heatmap

strat_plot_seg <- file.path(cnv_plot_saving_dir, "cnv_seg_heatmap_split")
one_off_plot_seg <- file.path(cnv_plot_saving_dir, "cnv_seg_heatmap_all")


h1 <- plot_heatmap(tf_list_seg_wide, 1,
  annotation_name_side = "left",
  legend_title = "z-score\n(Log2Ratio, seg)",
  legend_at = seq(-0.5, 0.5, 0.25),
  show_detailed_sample_annotation = FALSE
)

h2 <- plot_heatmap(tf_list_seg_wide, 2,
  show_heatmap_legend = FALSE,
  legend_at = seq(-0.5, 0.5, 0.25),
  show_chr_name = FALSE,
  show_detailed_sample_annotation = FALSE,
  show_annotation_name = FALSE,
  show_point_anno_axis = FALSE
)

h3 <- plot_heatmap(tf_list_seg_wide, 3,
  show_heatmap_legend = FALSE,
  legend_at = seq(-0.5, 0.5, 0.25),
  show_chr_name = FALSE,
  show_detailed_sample_annotation = FALSE,
  show_annotation_name = FALSE,
  show_point_anno_axis = FALSE
)

h4 <- plot_heatmap(tf_list_seg_wide, 4,
  show_heatmap_legend = FALSE,
  legend_at = seq(-0.5, 0.5, 0.25),
  show_chr_name = FALSE,
  show_detailed_sample_annotation = FALSE,
  show_annotation_name = FALSE,
  show_point_anno_axis = FALSE
)

# h5 <- plot_heatmap(tf_list_seg_wide, 5,
#   show_heatmap_legend = FALSE,
#   legend_at = seq(-0.5, 0.5, 0.25),
#   show_chr_name = FALSE,
#   show_detailed_sample_annotation = FALSE,
#   show_annotation_name = FALSE,
#   show_point_anno_axis = FALSE
# )

h_combine_seg <- h1 + h2 + h3 + h4
# plot as heatmap


pdf(paste0(strat_plot_seg, ".pdf"), width = 4.88, height = 3.78)
draw(h_combine_seg, ht_gap = unit(0.2, "cm"), merge_legend = TRUE)
dev.off()




# plot all samples at a time ----------------------------------------------


# combine everything
# bind the matrx by columns and convert to list
tf_list_seg_wide_all <- purrr::reduce(tf_list_seg_wide, cbind) %>% list()
hall_seg <- plot_heatmap(tf_list_seg_wide_all, 1,
  annotation_name_side = "left",
  legend_title = "z-score\n(Log2Ratio, seg)",
  legend_at = seq(-0.5, 0.5, 0.25)
)


# pdf(paste0(one_off_plot, ".pdf"), width = 4.88, height = 3.78)
pdf(paste0(one_off_plot_seg, ".pdf"), width = 18, height = 9.6)
draw(hall_seg)
dev.off()


png(paste0(one_off_plot_seg, ".png"), width = 18, height = 9.6, units = "in", res = 300)
draw(hall_seg)
dev.off()

################################################################################
# plot median sd for each bicohor
################################################################################

all_cnv_med <- all_cnv %>%
  group_by(bicohort, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

p_cnv_med <- ggplot(all_cnv_med, aes(x, y, group = bicohort)) +
  geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.4) +
  ylab("Median CNV") +
  xlab("Bin size") +
  theme_classic()



ggsave(
  filename = file.path(cnv_plot_saving_dir, "all_cnv_med.pdf"),
  plot = p_cnv_med,
  width = 10,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(cnv_plot_saving_dir, "all_cnv_med.pdf"))

# plot cnv, facet by author
p_cnv <- ggplot(all_cnv, aes(x, y, group = primary)) +
  geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
  xlab("Bin (5MB)") +
  ylab("CNV") +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  theme_classic() +
  # remove legend title
  guides(color = guide_legend(title = NULL)) +
  # set color palette
  scale_color_nord(palette = "aurora") +
  facet_wrap(~author)

ggsave(
  filename = file.path(cnv_plot_saving_dir, "all_cnv_facet.pdf"),
  plot = p_cnv,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(cnv_plot_saving_dir, "all_cnv_facet.pdf"))



# plot sd, facet by author and bicohort
p_cnv <- ggplot(all_cnv, aes(x, y, group = primary)) +
  geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
  xlab("Bin (5MB)") +
  ylab("CNV") +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  theme_classic() +
  # hide legend
  # theme(legend.position = "none") +
  # remove legend title
  guides(color = guide_legend(title = NULL)) +
  # set color palette
  scale_color_manual(values = c(
    "Cancer" = "grey",
    "Healthy" = "royalblue"
  )) +
  facet_grid(bicohort ~ author)

ggsave(
  filename = file.path(cnv_plot_saving_dir, "all_cnv_facet_author_cohort.pdf"),
  plot = p_cnv,
  width = 10,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(cnv_plot_saving_dir, "all_cnv_facet_author_cohort.pdf"))


################################################################################
# plot cnv facet by ichorcna_tf_strat
################################################################################


# count the number of sample for each ichorcna_tf_strat
make_label <- all_cnv %>%
  group_by(ichorcna_tf_strat) %>%
  summarize(n = length(unique(primary))) %>%
  ungroup() %>%
  mutate(facet_label = paste0("ichorCNA TF ", ichorcna_tf_strat, "\n", "(n=", n, ")"))

facet_label <- make_label$facet_label
names(facet_label) <- make_label$ichorcna_tf_strat


((
  p_cnv <- ggplot(all_cnv, aes(x, y, group = primary)) +
    # geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = 1, fill = "lightgrey", alpha = 0.2)+
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.8, linewidth = 0.2) +
    xlab("Bin index (5 MB)") +
    ylab("CNV") +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # make x axis text smaller
    theme(axis.text.x = element_text(size = 2)) +
    # set color palette
    scale_color_nord(palette = "aurora") +
    facet_wrap(~ichorcna_tf_strat, labeller = labeller(ichorcna_tf_strat = facet_label)) +
    # set facet label backgroud color to light grey
    theme(strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"))

))

plot_file <- file.path(cnv_plot_saving_dir, "all_cnv_facet_ichorcna_tf_strat.pdf")
# save plot
ggsave(
  filename = plot_file,
  plot = p_cnv,
  width = 7,
  height = 4,
  dpi = 300
)
message("Saved to ", plot_file)


################################################################################
# also use plotly to save p_long as html file
################################################################################
p_html <- ggplotly(p_cnv)

# save to htm
htmlwidgets::saveWidget(p_html, file.path(cnv_plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.html"), selfcontained = TRUE)


################################################################################
# boxplot: correlation with median healthy (non-seg) for each ichorcna_tf_strat
################################################################################
median_healthy <- all_cnv %>%
  filter(cohort == "Healthy") %>%
  group_by(x) %>%
  summarize(y = median(y)) %>%
  mutate(chr = str_extract(x, "\\d+")) %>%
  mutate(chr = as.numeric(chr)) %>%
  mutate(arm = str_extract(x, "[pq]")) %>%
  mutate(arm_index = str_extract(x, "_(\\d+)", group = 1)) %>%
  mutate(arm_index = as.numeric(arm_index)) %>%
  arrange(chr, arm, arm_index)

median_healthy$x <- factor(median_healthy$x, levels = median_healthy$x)
all_cnv$x <- factor(all_cnv$x, levels = levels(median_healthy$x))


all_cnv_cor <- all_cnv %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarise(cor_with_median_healthy = cor(y, median_healthy$y))

median_h_value <- all_cnv_cor %>%
  filter(ichorcna_tf_strat == "Healthy") %>%
  pull(cor_with_median_healthy) %>%
  median()


# plot as box plot
ylab_txt <- paste("Cor with median healthy CNV", sep = "")
p_box_cnv_cor_with_healthy <- ggplot(all_cnv_cor, aes(x = ichorcna_tf_strat, y = cor_with_median_healthy)) +
  geom_violin(fill = "lightgrey", color = "lightgrey", width = 1) +
  # boxplot without outliers
  geom_jitter(
    color = "black",
    fill = "#e9e6e6",
    alpha = 0.9,
    width = 0.2,
    size = 0.2,
    # shape = 21,
    stroke = 0.01
  ) +
  geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA, width = 0.15, linewidth = 0.2) +
  ylab(ylab_txt) +
  xlab("ichorCNA Tumor Fraction") +
  geom_hline(yintercept = median_h_value, color = "black", linetype = "dashed", linewidth = 0.1) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.2)) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  theme(axis.title = element_text(size = 5.5), axis.text = element_text(size = 5)) +
  theme(axis.title.x = element_blank())

plot_file <- file.path(cnv_plot_saving_dir, "all_cnv_cor_with_healthy.pdf")

ggsave(
  filename = plot_file,
  plot = p_box_cnv_cor_with_healthy,
  width = 4.27,
  height = 4.72,
  units = "cm",
  dpi = 300
)

message("Saved to ", plot_file)
################################################################################
# boxplot: correlation with median healthy (segmented) for each ichorcna_tf_strat
################################################################################
message("start to plot correlation with median healthy (segmented)")
median_healthy <- all_cnv_seg %>%
  filter(cohort == "Healthy") %>%
  group_by(x) %>%
  summarize(y = median(y)) %>%
  mutate(chr = str_extract(x, "\\d+")) %>%
  mutate(chr = as.numeric(chr)) %>%
  mutate(arm = str_extract(x, "[pq]")) %>%
  mutate(arm_index = str_extract(x, "_(\\d+)", group = 1)) %>%
  mutate(arm_index = as.numeric(arm_index)) %>%
  arrange(chr, arm, arm_index)

median_healthy$x <- factor(median_healthy$x, levels = median_healthy$x)
all_cnv_seg$x <- factor(all_cnv_seg$x, levels = levels(median_healthy$x))


all_cnv_cor <- all_cnv_seg %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarise(cor_with_median_healthy = cor(y, median_healthy$y))

median_h_value <- all_cnv_cor %>%
  filter(ichorcna_tf_strat == "Healthy") %>%
  pull(cor_with_median_healthy) %>%
  median()


# plot as box plot
ylab_txt <- paste("Cor with median healthy CNV", sep = "")
p_box_cnv_cor_with_healthy <- ggplot(all_cnv_cor, aes(x = ichorcna_tf_strat, y = cor_with_median_healthy)) +
  geom_violin(fill = "lightgrey", color = "lightgrey", width = 1) +
  # boxplot without outliers
  geom_jitter(
    color = "black",
    fill = "#e9e6e6",
    alpha = 0.9,
    width = 0.2,
    size = 0.2,
    # shape = 21,
    stroke = 0.01
  ) +
  geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA, width = 0.15, linewidth = 0.2) +
  ylab(ylab_txt) +
  xlab("ichorCNA Tumor Fraction") +
  geom_hline(yintercept = median_h_value, color = "black", linetype = "dashed", linewidth = 0.1) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.2)) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  theme(axis.title = element_text(size = 5.5), axis.text = element_text(size = 5)) +
  theme(axis.title.x = element_blank())

plot_file <- file.path(cnv_plot_saving_dir, "all_cnv_cor_with_healthy_segmented.pdf")

ggsave(
  filename = plot_file,
  plot = p_box_cnv_cor_with_healthy,
  width = 4.27,
  height = 4.72,
  units = "cm",
  dpi = 300
)

message("Saved to ", plot_file)


# ################################################################################
# # barplot: sd of sd for each ichorcna_tf_strat
# ################################################################################

# sd_sd_func <- function(input, min, max) {
#   ans <- input %>%
#     filter(x >= !!min & x <= !!max) %>%
#     group_by(primary, ichorcna_tf_strat) %>%
#     summarize(y_sd = sd(y)) %>%
#     ungroup()
#   return(ans)
# }

# all_sd_healthy_sd_range <- sd_sd_func(all_sd_healthy, min, max)
# all_sd_cancer_sd_range <- sd_sd_func(all_sd_cancer, min, max)

# # combine the two dataframes
# all_sd_sd_range <- rbind(all_sd_healthy_sd_range, all_sd_cancer_sd_range)

# # set the order of ichorcna_tf_strat, "Healthy" first
# all_sd_sd_range$ichorcna_tf_strat <- factor(all_sd_sd_range$ichorcna_tf_strat,
#   levels = label_levels
# )

# # plot as box plot
# sd_sd_ylab_txt <- paste("SD of SD (", min, "~", max, "bp fragments)", sep = "")
# ((
#   p_box_sd_sd <- ggplot(all_sd_sd_range, aes(x = ichorcna_tf_strat, y = y_sd)) +
#     # boxplot without outliers
#     geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA) +
#     geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.55) +
#     ylab(sd_sd_ylab_txt) +
#     xlab("ichorCNA Tumor Fraction") +
#     stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy") +
#     theme_classic() +
#     # remove legend
#     theme(legend.position = "none") +
#     scale_fill_nord(palette = "aurora") +
#     # rotate x-axis tick and text by 45 degree
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ))


# ggsave(
#   filename = file.path(cnv_plot_saving_dir, paste("all_sd_", min, "_", max, "_sd_sd_boxplot.pdf", sep = "")),
#   plot = p_box_sd_sd,
#   width = 6,
#   height = 4,
#   dpi = 300
# )

# message("Saved to ", file.path(cnv_plot_saving_dir, paste("all_sd_", min, "_", max, "_sd_sd_boxplot.pdf", sep = "")))

# ################################################################################
# # median sd for each ichorcna_tf_strat
# ################################################################################

# # calculate the median of all samples
# all_sd_healthy_med <- all_sd_healthy %>%
#   group_by(x) %>%
#   summarize(y = median(y)) %>%
#   ungroup() %>%
#   mutate(ichorcna_tf_strat = "Healthy")


# all_sd_cancer_med <- all_sd_cancer %>%
#   group_by(ichorcna_tf_strat, x) %>%
#   summarize(y = median(y)) %>%
#   ungroup()

# all_sd_med <- rbind(all_sd_healthy_med, all_sd_cancer_med)

# # set the order of ichorcna_tf_strat, "Healthy" first
# all_sd_med$ichorcna_tf_strat <- factor(all_sd_med$ichorcna_tf_strat,
#   levels = label_levels
# )

# # plot median sd for each ichorcna_tf_strat
# ((
#   p_long_med <- ggplot(all_sd_med, aes(x, y, group = ichorcna_tf_strat)) +
#     geom_rect(xmin = min, xmax = max, ymin = 0, ymax = 1, fill = "lightgrey", alpha = 0.2) +
#     geom_line(aes(color = ichorcna_tf_strat), alpha = 0.95, linewidth = 0.7) +
#     geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
#     ylab("Median SD") +
#     xlab("Fragment Length (bp)") +
#     # add 155 and 167 bp to the x-axis tick and text
#     scale_x_continuous(breaks = c(50, 80, 100, 155, 167, 200, 250)) +
#     # only show x between 80 and 250
#     coord_cartesian(xlim = c(80, 250)) +
#     # remove legend title
#     guides(color = guide_legend(title = NULL)) +
#     # increase the legend text size and line size
#     theme(
#       legend.text = element_text(size = 8),
#       legend.title = element_text(size = 8),
#       legend.key.height = unit(0.9, "cm")
#     ) +
#     theme_classic() +
#     scale_color_nord(palette = "aurora") +
#     # rotate x-axis tick and text by 45 degree
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ))



# ggsave(
#   filename = file.path(cnv_plot_saving_dir, "all_sd_med_ichorcna_tf_strat.pdf"),
#   plot = p_long_med,
#   width = 6,
#   height = 3,
#   dpi = 300
# )

# message("Saved to ", file.path(cnv_plot_saving_dir, "all_sd_med_ichorcna_tf_strat.pdf"))


# ################################################################################
# # use patchwork to combine the  plots, hide all legends
# ################################################################################

# p_long_med_nolegend <- p_long_med + theme(legend.position = "none")

# ctdna_is_shorter <- p_long +
#   (p_long_med_nolegend / p_box_sum_sd) +
#   plot_layout(widths = c(1.3, 1))

# ggsave(
#   filename = file.path(cnv_plot_saving_dir, "ctdna_has_higher_sd.pdf"),
#   plot = ctdna_is_shorter,
#   width = 10,
#   height = 6,
#   dpi = 300
# )
# message("Saved to ", file.path(cnv_plot_saving_dir, "ctdna_has_higher_sd.pdf"))
