library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(ComplexHeatmap)

sl_plot_saving_dir_base <- "/home/nrlab/wang04/ulyses/plot_mae/plots/sl_ratio"

source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")

sl_plot_saving_dir <- file.path(sl_plot_saving_dir_base, paste0("iteration_", iteration_num_hyper))

# if the directory does not exist, create it
if (!dir.exists(sl_plot_saving_dir)) {
  dir.create(sl_plot_saving_dir)
}


################################################################################
# sl_ratio
################################################################################
all_sl <- all_long_filter %>%
  dplyr::filter(assay %in% c("sl_ratio")) %>%
  # remove rows 6p_6
  dplyr::filter(rowname != "6p_6") %>%
  dplyr::mutate(chr_arm = str_extract(rowname, "\\d+[pq]")) %>%
  dplyr::mutate(chr = str_extract(chr_arm, "\\d+")) %>%
  dplyr::mutate(arm = str_extract(chr_arm, "[pq]")) %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)

all_sl_cnv <- all_long_filter %>%
  dplyr::filter(assay %in% c("sl_ratio", "cnv")) %>%
  # remove rows 6p_6
  dplyr::filter(rowname != "6p_6") %>%
  dplyr::mutate(chr_arm = str_extract(rowname, "\\d+[pq]")) %>%
  dplyr::mutate(chr = str_extract(chr_arm, "\\d+")) %>%
  dplyr::mutate(arm = str_extract(chr_arm, "[pq]")) %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)



# plot median sd for each bicohor

all_sl_med <- all_sl %>%
  group_by(bicohort, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

((
  p_sl_med <- ggplot(all_sl_med, aes(x, y, group = bicohort)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.4) +
    ylab("Median S/L Ratio") +
    xlab("Bin size") +
    scale_color_manual(values = c(
      "Healthy" = healthy_color_hyper,
      "Cancer" = cancer_color_hyper
    )) +
    theme_classic()

))


ggsave(
  filename = file.path(sl_plot_saving_dir, "all_sl_med.pdf"),
  plot = p_sl_med,
  width = 10,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(sl_plot_saving_dir, "all_sl_med.pdf"))

# plot s/l, facet by author
((
  p_sl <- ggplot(all_sl, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
    xlab("Bin (5MB)") +
    ylab("S/L Ratio") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
    theme_classic() +
    # remove legend title
    guides(color = guide_legend(title = NULL)) +
    scale_color_manual(values = c(
      "Healthy" = healthy_color_hyper,
      "Cancer" = cancer_color_hyper
    )) +
    facet_wrap(~author)
))

ggsave(
  filename = file.path(sl_plot_saving_dir, "all_sl_facet.pdf"),
  plot = p_sl,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(sl_plot_saving_dir, "all_sl_facet.pdf"))



# plot sd, facet by author and bicohort
((
  p_sl <- ggplot(all_sl, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
    xlab("Bin (5MB)") +
    ylab("S/L Ratio") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
    theme_classic() +
    # hide legend
    # theme(legend.position = "none") +
    # remove legend title
    guides(color = guide_legend(title = NULL)) +
    # set color palette
    scale_color_manual(values = c(
      "Healthy" = healthy_color_hyper,
      "Cancer" = cancer_color_hyper
    )) +
    facet_grid(bicohort ~ author)
))

ggsave(
  filename = file.path(sl_plot_saving_dir, "all_sl_facet_author_cohort.pdf"),
  plot = p_sl,
  width = 10,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(sl_plot_saving_dir, "all_sl_facet_author_cohort.pdf"))


################################################################################
# plot sl facet by ichorcna_tf_strat
################################################################################


# count the number of sample for each ichorcna_tf_strat
make_label <- all_sl %>%
  group_by(ichorcna_tf_strat) %>%
  summarize(n = length(unique(primary))) %>%
  ungroup() %>%
  mutate(facet_label = paste0("ichorCNA TF ", ichorcna_tf_strat, "\n", "(n=", n, ")"))

facet_label <- make_label$facet_label
names(facet_label) <- make_label$ichorcna_tf_strat


((
  p_sl <- ggplot(all_sl, aes(x, y, group = primary)) +
    # geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = 1, fill = "lightgrey", alpha = 0.2)+
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.8, linewidth = 0.2) +
    xlab("Bin index (5 MB)") +
    ylab("S/L Ratio") +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # make x axis text smaller
    theme(axis.text.x = element_text(size = 2)) +
    # set color palette
    scale_color_manual(values = ichorcna_tf_strat_color_hyper[unique(all_sl$ichorcna_tf_strat) |> sort()]) +
    facet_wrap(~ichorcna_tf_strat, labeller = labeller(ichorcna_tf_strat = facet_label)) +
    # set facet label backgroud color to light grey
    theme(strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"))

))

plot_file <- file.path(sl_plot_saving_dir, "all_sl_facet_ichorcna_tf_strat.pdf")
# save plot
ggsave(
  filename = plot_file,
  plot = p_sl,
  width = 7,
  height = 4,
  dpi = 300
)
message("Saved to ", plot_file)


# also use plotly to save p_long as html file
p_html <- ggplotly(p_sl)

# save to htm
htmlwidgets::saveWidget(p_html, file.path(sl_plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.html"), selfcontained = TRUE)


################################################################################
# boxplot: correlation with median healthy for each ichorcna_tf_strat
################################################################################
median_healthy <- all_sl %>%
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
all_sl$x <- factor(all_sl$x, levels = levels(median_healthy$x))


all_sl_cor <- all_sl %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarise(cor_with_median_healthy = cor(y, median_healthy$y))

# calcualte the median healthy S/L ratio
median_h_value <- all_sl_cor %>%
  filter(ichorcna_tf_strat == "Healthy") %>%
  pull(cor_with_median_healthy) %>%
  median()
# plot as box plot
ylab_txt <- paste("Cor with median healthy S/L Ratio", sep = "")
p_box_sl_cor_with_healthy <- ggplot(all_sl_cor, aes(x = ichorcna_tf_strat, y = cor_with_median_healthy)) +
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
  # boxplot without outliers
  geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA, width = 0.15, linewidth = 0.2) +
  geom_hline(yintercept = median_h_value, color = "black", linetype = "dashed", linewidth = 0.1) +
  ylab(ylab_txt) +
  xlab("ichorCNA Tumor Fraction") +
  # add p-value to the plot
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  theme(axis.title = element_text(size = 5.5), axis.text = element_text(size = 5)) +
  theme(axis.title.x = element_blank()) 


plot_file <- file.path(sl_plot_saving_dir, "all_sl_cor_with_healthy.pdf")

ggsave(
  filename = plot_file,
  plot = p_box_sl_cor_with_healthy,
  width = 4.27,
  height = 4.72,
  units = "cm",
  dpi = 300
)

message("Saved to ", plot_file)


################################################################################
# Nitzan: barplot: correlation with median Breast for each ichorcna_tf_strat
################################################################################
median_breast <- all_sl %>%
  filter(cohort == "Breast") %>%
  group_by(x) %>%
  summarize(y = median(y)) %>%
  mutate(chr = str_extract(x, "\\d+")) %>%
  mutate(chr = as.numeric(chr)) %>%
  mutate(arm = str_extract(x, "[pq]")) %>%
  mutate(arm_index = str_extract(x, "_(\\d+)", group = 1)) %>%
  mutate(arm_index = as.numeric(arm_index)) %>%
  arrange(chr, arm, arm_index)

median_breast$x <- factor(median_breast$x, levels = median_breast$x)
all_sl$x <- factor(all_sl$x, levels = levels(median_breast$x))


all_sl_cor <- all_sl %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarise(cor_with_median_breast = cor(y, median_breast$y))

all_sl_cor <- all_sl %>%
  group_by(primary, ichorcna_tf_strat, cohort) %>%
  summarise(cor_with_median_breast = cor(y, median_breast$y))

# plot as box plot
ylab_txt <- paste("Cor with median Breast S/L Ratio", sep = "")
((
  p_box_sl_cor_with_breast <- ggplot(all_sl_cor, aes(x = cohort, y = cor_with_median_breast)) +
    # boxplot without outliers
    geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA) +
    geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.55) +
    ylab(ylab_txt) +
    xlab("ichorCNA Tumor Fraction") +
    # add p-value to the plot
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Breast") +
    theme_classic() +
    # remove legend
    theme(legend.position = "none") +
    scale_fill_manual(values = ichorcna_tf_strat_color_hyper[unique(all_sl_cor$ichorcna_tf_strat) |> sort()]) +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ichorcna_tf_strat)

))

plot_file <- file.path(sl_plot_saving_dir, "all_sl_cor_with_breast_facet_by_ichorcna_tf_strat.pdf")

ggsave(
  filename = plot_file,
  plot = p_box_sl_cor_with_breast,
  width = 18,
  height = 7,
  dpi = 300
)

message("Saved to ", plot_file)





################################################################################
# plot sl_ratio as heatmap using ComplexHeatmap
################################################################################
# plot as heatmap
# use ComplexHeatmap package
# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html




all_sl_grouped <- all_sl %>%
  group_by(ichorcna_tf_strat)

group_keys(all_sl_grouped)

tf_list <- all_sl %>%
  group_by(ichorcna_tf_strat) %>%
  group_split()

names(tf_list) <- group_keys(all_sl_grouped) |> pull(ichorcna_tf_strat)

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

saveRDS(tf_list_wide, file.path(sl_plot_saving_dir, "tf_list_wide.rds"))

# plot as heatmap


library(ComplexHeatmap)

plot_heatmap <- function(input, index, show_heatmap_legend = TRUE) {
  row_title <- paste( "5 Mb Genomic Bins ", "(n=", dim(input[[index]])[[1]], ")", sep = "")

  col_title <- paste(names(input)[[index]], "\n (", "n=", dim(input[[index]])[[2]], ")")
  p <- Heatmap(tf_list_wide[[index]],
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_title = row_title,
    row_title_gp = gpar(fontsize = 5),
    column_title = col_title,
    column_title_gp = gpar(fontsize = 5),
    show_heatmap_legend = show_heatmap_legend,
    heatmap_legend_param = list(title = "z-score", 
                                at = seq(-5, 5, 2),
                                legend_height = unit(2, "cm"),
                                grid_width = unit(0.13, "cm"),
                                # set title size to 5
                                title_gp = gpar(fontsize = 5),
                                labels_gp = gpar(fontsize = 5)),
    use_raster = FALSE,
    raster_device = c("png"),
    raster_quality = 30,
  )

  return(p)
}


h1 <- plot_heatmap(tf_list_wide, 1)
h2 <- plot_heatmap(tf_list_wide, 2, show_heatmap_legend = FALSE)
h3 <- plot_heatmap(tf_list_wide, 3, show_heatmap_legend = FALSE)
h4 <- plot_heatmap(tf_list_wide, 4, show_heatmap_legend = FALSE)
# h5 <- plot_heatmap(tf_list_wide, 5, show_heatmap_legend = FALSE)

h_combine <- h1 + h2 + h3 + h4



sl_heatmap_file <- file.path(sl_plot_saving_dir, "sl_heatmap_all.png")
png(sl_heatmap_file, width =83.5, height = 44.8, units = "mm", res = 300)
draw(h_combine, ht_gap = unit(0.06, "cm"), merge_legend = TRUE)
dev.off()

pdf(paste0(sl_heatmap_file, ".pdf"), width = 3.287402, height = 1.748031)
draw(h_combine, ht_gap = unit(0.06, "cm"), merge_legend = TRUE)
dev.off()
