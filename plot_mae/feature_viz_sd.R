library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(factoextra)


plot_saving_dir_base <- "/home/nrlab/wang04/ulyses/plot_mae/plots/sd"
source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")

# make plot_saving_dir with iteration_num_hyper
plot_saving_dir <- file.path(plot_saving_dir_base, paste0("iteration_", iteration_num_hyper))
if (!dir.exists(plot_saving_dir)) {
  dir.create(plot_saving_dir)
}


################################################################################
# length
################################################################################
all_sd <- all_long_filter %>%
  dplyr::filter(assay == "sd") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)


# plot all samples together, color by author
((
  p_long <- ggplot(all_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.4, linewidth = 0.25) +
    ylab("SD") +
    xlab("Fragment Length (bp)") +
    theme_classic()

))

ggsave(
  filename = file.path(plot_saving_dir, "all_sd.pdf"),
  plot = p_long,
  width = 16,
  height = 7,
  dpi = 300,
  units = "cm"
)
message("Saved to ", file.path(plot_saving_dir, "all_sd.pdf"))

# plot healthy and cancer seperately


data_selected <- all_sd %>% dplyr::filter(bicohort %in% c("Healthy"))

n_sample <- data_selected |>
  pluck("primary") |>
  unique() |>
  length()



p_healthy <- ggplot(data_selected, aes(x, y, group = primary)) +
  geom_line(color = "grey", alpha = 0.4, linewidth = 0.25) +
  ylab("SD") +
  xlab("Fragment Length (bp)") +
  # add n_sample as the title
  ggtitle(paste0("Healthy (n=", n_sample, ")")) +
  # set y between 0 and 0.05
  # coord_cartesian(ylim = c(0, 0.05)) +
  # add vertical dashed lines at 167 bp
  # geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  # add mean and max x to x-axis tick and text
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  theme_classic()


data_selected <- all_sd %>% dplyr::filter(bicohort %in% c("Cancer"))

n_sample <- data_selected |>
  pluck("primary") |>
  unique() |>
  length()

p_cancer <- ggplot(data_selected, aes(x, y, group = primary)) +
  geom_line(color = "red", alpha = 0.4, linewidth = 0.25) +
  ylab("SD") +
  xlab("Fragment Length (bp)") +
  # add n_sample as the title
  ggtitle(paste0("Cancer (n=", n_sample, ")")) +
  # set y between 0 and 0.05
  # coord_cartesian(ylim = c(0, 0.05)) +
  # add vertical dashed lines at 167 bp
  # geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  # add mean and max x to x-axis tick and text
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  theme_classic()

p_healthy_cancer_sep <- p_healthy / p_cancer



ggsave(
  filename = file.path(plot_saving_dir, "all_sd_cancer_healthy_stacked.pdf"),
  plot = p_healthy_cancer_sep,
  width = 16,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", file.path(plot_saving_dir, "all_sd_cancer_healthy_stacked.pdf"))


# plot without primary starting with "EE"

data_selected <- all_sd %>%
  dplyr::filter(bicohort == "Healthy") %>%
  dplyr::filter(!str_detect(primary, "^EE"))

n_sample <- data_selected |>
  pluck("primary") |>
  unique() |>
  length()


p_healthy <- ggplot(data_selected, aes(x, y, group = primary)) +
  geom_line(color = "grey", alpha = 0.4, linewidth = 0.25) +
  ylab("SD") +
  xlab("Fragment Length (bp)") +
  # add n_sample as the title
  ggtitle(paste0("Healthy excluding FinaleDB (n=", n_sample, ")")) +
  # set y between 0 and 0.05
  # coord_cartesian(ylim = c(0, 0.05)) +
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  theme_classic()

data_selected <- all_sd %>%
  dplyr::filter(bicohort == "Cancer") %>%
  dplyr::filter(!str_detect(primary, "^EE"))

n_sample <- data_selected |>
  pluck("primary") |>
  unique() |>
  length()

p_cancer <- ggplot(data_selected, aes(x, y, group = primary)) +
  geom_line(color = "red", alpha = 0.4, linewidth = 0.25) +
  ylab("SD") +
  xlab("Fragment Length (bp)") +
  # add n_sample as the title
  ggtitle(paste0("Cancer excluding FinaleDB (n=", n_sample, ")")) +
  # set y between 0 and 0.05
  # coord_cartesian(ylim = c(0, 0.05)) +
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  theme_classic()

p_healthy_cancer_sep <- p_healthy / p_cancer



ggsave(
  filename = file.path(plot_saving_dir, "all_sd_cancer_healthy_stacked_without_FINALEDB.pdf"),
  plot = p_healthy_cancer_sep,
  width = 16,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", file.path(plot_saving_dir, "all_sd_cancer_healthy_stacked_without_FINALEDB.pdf"))

p_healthy_html <- ggplotly(p_healthy)

# save to html
htmlwidgets::saveWidget(p_healthy_html, file.path(plot_saving_dir, "all_sd_healthy_without_FINALEDB.html"))

message("Saved to", file.path(plot_saving_dir, "all_sd_healthy_without_FINALEDB.html"))




# plot median length for each bicohor

all_sd_med <- all_sd %>%
  group_by(bicohort, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

((
  p_long_med <- ggplot(all_sd_med, aes(x, y, group = bicohort)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.4) +
    ylab("Median SD") +
    xlab("Fragment Length (bp)") +
    scale_x_continuous(breaks = c(min(all_sd_med$x), 100, 167, 200, max(all_sd_med$x))) +
    theme_classic()

))


ggsave(
  filename = file.path(plot_saving_dir, "all_sd_med.pdf"),
  plot = p_long_med,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_sd_med.pdf"))






###############################################################################
# do PCA analysis
###############################################################################

###############################################################################
# pca ONLY healthy
###############################################################################
data_selected_sd_cancer <- all_sd %>% dplyr::filter(bicohort %in% c("Cancer"))
data_selected_sd_healthy <- all_sd %>% dplyr::filter(bicohort %in% c("Healthy"))

all_sd1 <- select(data_selected_sd_healthy, primary, rowname, value, author, data_source, bicohort)
# expand the "rowname" columns to wider format
all_sd2 <- pivot_wider(all_sd1, names_from = rowname, values_from = value) %>%
  # change primary to rownames
  column_to_rownames("primary")
all_sd2_active <- all_sd2[, 4:ncol(all_sd2)]
res.pca <- prcomp(all_sd2_active)


length_range_label_isize_min <- min(all_sd1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label_isize_max <- max(all_sd1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label <- paste0("(", length_range_label_isize_min, " - ", length_range_label_isize_max, "bp)")
pca_plot_title <- paste0("PCA of Fragment Length ", length_range_label)

# plot file suffix
plot_suffix <- "_pca_sd_healthy.pdf"

fviz_eig(res.pca)
groups_bicohort <- all_sd2$bicohort
groups_bicohort <- factor(groups_bicohort, levels = c("Healthy", "Cancer"))
groups_author <- all_sd2$author
groups_data_source <- all_sd2$data_source



p_pca_sd_bicohort <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_bicohort, # color by groups
  palette = c("#00AFBB", "red"),
  alpha.ind = 0.45,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  legend.title = "Cohort",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_bicohort", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_bicohort,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)




#------------------------------------------------------------------------------

authors_vec <- all_sd2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_sd_author <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_author, # color by groups
  palette = authors_color_code,
  alpha.ind = 0.45,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  title = pca_plot_title,
  legend.title = "Author",
  repel = TRUE
)
p_pca_sd_author <- p_pca_sd_author + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(plot_saving_dir, paste0("all_sd_author", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_author,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

#------------------------------------------------------------------------------

p_pca_sd_datasource <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_data_source, # color by groups
  palette = c("grey", "royalblue", "orange3"),
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Data Source",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_datasource", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_datasource,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)


###############################################################################
# pca: healthy + cancer
###############################################################################
data_selected_sd_cancer <- all_sd %>% dplyr::filter(bicohort %in% c("Cancer"))
data_selected_sd_healthy <- all_sd %>% dplyr::filter(bicohort %in% c("Healthy"))

all_sd1 <- select(all_sd, primary, rowname, value, author, data_source, bicohort, ichorcna_tf_strat)
# expand the "rowname" columns to wider format
all_sd2 <- pivot_wider(all_sd1, names_from = rowname, values_from = value) %>%
  # change primary to rownames
  column_to_rownames("primary")
all_sd2_active <- all_sd2[, 5:ncol(all_sd2)]
res.pca <- prcomp(all_sd2_active)


length_range_label_isize_min <- min(all_sd1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label_isize_max <- max(all_sd1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label <- paste0("(", length_range_label_isize_min, " - ", length_range_label_isize_max, "bp)")
pca_plot_title <- paste0("PCA of Fragment Length ", length_range_label)

# plot file suffix
plot_suffix <- "_pca_sd_healthy_and_cancer.pdf"

# fviz_eig(res.pca)
groups_bicohort <- all_sd2$bicohort
groups_bicohort <- factor(groups_bicohort, levels = c("Healthy", "Cancer"))
groups_author <- all_sd2$author
groups_data_source <- all_sd2$data_source
groups_ichorcna_tf_strat <- all_sd2$ichorcna_tf_strat

#--------------------------------------------------------------------------
# plot by ichorcna_tf_strat
p_pca_sd_ichorcna_tf_strat <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_ichorcna_tf_strat, # color by groups
  palette = ichorcna_tf_strat_color_hyper[groups_ichorcna_tf_strat],
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "ichorCNA TF Strat",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_ichorcna_tf_strat", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_ichorcna_tf_strat,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)



p_pca_sd_bicohort <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_bicohort, # color by groups
  # palette = c("royalblue", "#f95959"),
  palette = c(healthy_color_hyper, cancer_color_hyper),
  alpha.ind = 0.4,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  legend.title = "Cohort",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_bicohort", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_bicohort,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)




#------------------------------------------------------------------------------

authors_vec <- all_sd2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_sd_author <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_author, # color by groups
  palette = authors_color_code,
  alpha.ind = 0.45,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  title = pca_plot_title,
  legend.title = "Author",
  repel = TRUE
)
p_pca_sd_author <- p_pca_sd_author + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(plot_saving_dir, paste0("all_sd_author", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_author,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

#------------------------------------------------------------------------------

p_pca_sd_datasource <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_data_source, # color by groups
  palette = c("grey", "royalblue", "orange3"),
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Data Source",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_datasource", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_datasource,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)



###############################################################################
# pca: cancer
###############################################################################
data_selected_sd_cancer <- all_sd %>% dplyr::filter(bicohort %in% c("Cancer"))
data_selected_sd_healthy <- all_sd %>% dplyr::filter(bicohort %in% c("Healthy"))

all_sd1 <- select(data_selected_sd_cancer, primary, rowname, value, author, data_source, bicohort, ichorcna_tf_strat)
# expand the "rowname" columns to wider format
all_sd2 <- pivot_wider(all_sd1, names_from = rowname, values_from = value) %>%
  # change primary to rownames
  column_to_rownames("primary")
all_sd2_active <- all_sd2[, 5:ncol(all_sd2)]
res.pca <- prcomp(all_sd2_active)


length_range_label_isize_min <- min(all_sd1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label_isize_max <- max(all_sd1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label <- paste0("(", length_range_label_isize_min, " - ", length_range_label_isize_max, "bp)")
pca_plot_title <- paste0("PCA of Fragment Length ", length_range_label)

# plot file suffix
plot_suffix <- "_pca_sd_cancer.pdf"

fviz_eig(res.pca)
groups_bicohort <- all_sd2$bicohort
groups_bicohort <- factor(groups_bicohort, levels = c("Healthy", "Cancer"))
groups_author <- all_sd2$author
groups_data_source <- all_sd2$data_source
groups_ichorcna_tf_strat <- all_sd2$ichorcna_tf_strat

#--------------------------------------------------------------------------
# plot by ichorcna_tf_strat
p_pca_sd_ichorcna_tf_strat <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_ichorcna_tf_strat, # color by groups
  palette = ichorcna_tf_strat_color_hyper[groups_ichorcna_tf_strat],
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "ichorCNA TF Strat",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_ichorcna_tf_strat", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_ichorcna_tf_strat,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)


#--------------------------------------------------------------------------
p_pca_sd_bicohort <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_bicohort, # color by groups
  palette = c("royalblue", "#f95959"),
  alpha.ind = 0.4,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  legend.title = "Cohort",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_bicohort", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_bicohort,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)




#------------------------------------------------------------------------------

authors_vec <- all_sd2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_sd_author <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_author, # color by groups
  palette = authors_color_code,
  alpha.ind = 0.45,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  title = pca_plot_title,
  legend.title = "Author",
  repel = TRUE
)
p_pca_sd_author <- p_pca_sd_author + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(plot_saving_dir, paste0("all_sd_author", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_author,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

#------------------------------------------------------------------------------

p_pca_sd_datasource <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_data_source, # color by groups
  palette = c("grey", "royalblue", "orange3"),
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Data Source",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(plot_saving_dir, paste0("all_sd_datasource", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_sd_datasource,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

################################################################################
# END of pca analysis
################################################################################

################################################################################
# plot length, facet by author

authors_vec <- all_sd$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_long <- ggplot(all_sd, aes(x, y, group = primary)) +
  geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
  xlab("Fragment Length (bp)") +
  ylab("SD") +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  theme_classic() +
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  # hide legend
  theme(legend.position = "none") +
  # set color palette
  # scale_color_nord(palette = "aurora") +
  scale_color_manual(values = c(
    "Cancer" = cancer_color_hyper,
    "Healthy" = healthy_color_hyper
  )) +
  # rotate x axis tick and text by 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~author) +
  # remove boarder of facet panels and set the background color to manual color values
  # remove boarder of facet panels
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "lightgrey"))


ggsave(
  filename = file.path(plot_saving_dir, "all_sd_facet.pdf"),
  plot = p_long,
  width = 12,
  height = 8,
  units = "cm",
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_sd_facet.pdf"))


# plot length, facet by author and bicohort
all_sd <- all_sd %>%
  mutate(bicohort = factor(bicohort, levels = c("Healthy", "Cancer")))

p_long <- ggplot(all_sd %>% filter(author %in% authors_vec[1:6]), aes(x, y, group = primary)) +
  geom_line(aes(color = bicohort), alpha = 0.45, linewidth = 0.18) +
  xlab("Fragment Length (bp)") +
  ylab("SD") +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  theme_classic() +
  # hide legend
  theme(legend.position = "none") +
  # set color palette
  scale_color_manual(values = c(
    "Cancer" = cancer_color_hyper,
    "Healthy" = healthy_color_hyper
  )) +
  # facet_grid(author ~ bicohort) +
  facet_grid(bicohort ~ author) +
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "lightgrey")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_long2 <- ggplot(all_sd %>% filter(author %in% authors_vec[7:length(authors_hyper)]), aes(x, y, group = primary)) +
  geom_line(aes(color = bicohort), alpha = 0.45, linewidth = 0.18) +
  xlab("Fragment Length (bp)") +
  ylab("SD") +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  scale_x_continuous(breaks = c(min(all_sd$x), 100, 167, 200, max(all_sd$x))) +
  theme_classic() +
  # hide legend
  theme(legend.position = "none") +
  # set color palette
  scale_color_manual(values = c(
    "Cancer" = cancer_color_hyper,
    "Healthy" = healthy_color_hyper
  )) +
  # facet_grid(author ~ bicohort) +
  facet_grid(bicohort ~ author) +
  theme(panel.border = element_blank()) +
  theme(strip.background = element_rect(fill = "lightgrey")) +
  # rotate x axis tick and text by 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plong_1_2 <- p_long / p_long2

ggsave(
  filename = file.path(plot_saving_dir, "all_sd_facet_author_cohort.pdf"),
  plot = plong_1_2,
  width = 18,
  height = 16,
  dpi = 300,
  units = "cm"
)

message("Saved to ", file.path(plot_saving_dir, "all_sd_facet_author_cohort.pdf"))


################################################################################
# main figure plot
################################################################################

# plot length facet by ichorcna_tf_strat
all_sd_healthy <- all_sd %>%
  filter(bicohort == "Healthy") %>%
  filter(ichorcna_tf_strat != "unknown") %>%
  mutate(ichorcna_tf_strat = "Healthy")

all_sd_cancer <- all_sd %>%
  filter(bicohort == "Cancer") %>%
  filter(ichorcna_tf_strat != "unknown")

# combine the two dataframes
all_sd_clean <- rbind(all_sd_healthy, all_sd_cancer)


label_levels <- ichorcna_tf_strat_hyper

# set the order of ichorcna_tf_strat, "Healthy" first
all_sd_clean$ichorcna_tf_strat <- factor(all_sd_clean$ichorcna_tf_strat,
  levels = label_levels
)

# count the number of sample for each ichorcna_tf_strat
make_label <- all_sd_clean %>%
  group_by(ichorcna_tf_strat) %>%
  summarize(n = length(unique(primary))) %>%
  ungroup() %>%
  mutate(facet_label = ifelse(ichorcna_tf_strat == "Healthy", paste0("Healthy", "\n", "(n=", n, ")"), paste0("ichorCNA TF ", ichorcna_tf_strat, "\n", "(n=", n, ")")))

facet_label <- make_label$facet_label
names(facet_label) <- make_label$ichorcna_tf_strat


p_long <- ggplot(all_sd_clean, aes(x, y, group = primary)) +
  # geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = Inf, fill = "lightgrey", alpha = 0.2)+
  geom_line(aes(color = ichorcna_tf_strat), alpha = 0.3, linewidth = 0.2) +
  xlab("Fragment Length (bp)") +
  ylab("SD") +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  theme_classic() +
  # hide legend
  theme(legend.position = "none") +
  # rotate x-axis tick and text by 45 degree
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_continuous(breaks = c(min(all_sd_clean$x), 100, 167, 200, max(all_sd_clean$x))) +
  # set line color manually to ichorcna_tf_strat_color_hyper[label_levels]
  scale_color_manual(values = ichorcna_tf_strat_color_hyper[label_levels]) +
  facet_wrap(~ichorcna_tf_strat, labeller = labeller(ichorcna_tf_strat = facet_label), nrow = 1) +
  # set facet label backgroud color to light grey
  theme(strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"))



# save plot
ggsave(
  filename = file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.pdf"),
  plot = p_long,
  width = 16,
  height = 9,
  dpi = 300,
  units = "cm"
)
message("Saved to ", file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.pdf"))


################################################################################
# also use plotly to save p_long as html file
################################################################################
p_html <- ggplotly(p_long)

# save to htm
htmlwidgets::saveWidget(p_html, file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.html"))
message("Saved to ", file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.html"))



# median length for each ichorcna_tf_strat

# calculate the median of all samples
all_sd_healthy_med <- all_sd_healthy %>%
  group_by(x) %>%
  summarize(y = median(y)) %>%
  ungroup() %>%
  mutate(ichorcna_tf_strat = "Healthy")


all_sd_cancer_med <- all_sd_cancer %>%
  group_by(ichorcna_tf_strat, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

all_sd_med <- rbind(all_sd_healthy_med, all_sd_cancer_med)

# set the order of ichorcna_tf_strat, "Healthy" first
all_sd_med$ichorcna_tf_strat <- factor(all_sd_med$ichorcna_tf_strat,
  levels = ichorcna_tf_strat_hyper
)

# plot median length for each ichorcna_tf_strat
p_long_med <- ggplot(all_sd_med, aes(x, y, group = ichorcna_tf_strat)) +
  geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = Inf, fill = "lightgrey", alpha = 0.2) +
  geom_line(aes(color = ichorcna_tf_strat), alpha = 0.95, linewidth = 0.4) +
  geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1) +
  ylab("Median SD") +
  xlab("Fragment Length (bp)") +
  # add 155 and 167 bp to the x-axis tick and text
  # scale_x_continuous(breaks = c(50, 80, 100, 155, 167, 200, 250)) +
  scale_x_continuous(breaks = c(min(all_sd_clean$x), 100, 155, 167, 200, max(all_sd_clean$x))) +
  # only show x between 80 and 250
  coord_cartesian(xlim = c(80, 250)) +
  # remove legend title
  guides(color = guide_legend(title = NULL)) +
  # increase the legend text size and line size
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.height = unit(0.9, "cm")
  ) +
  theme_classic() +
  # scale_color_nord(palette = "aurora") +
  scale_color_manual(values = ichorcna_tf_strat_color_hyper[label_levels]) +
  # rotate x-axis tick and text by 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggsave(
  filename = file.path(plot_saving_dir, "all_sd_med_ichorcna_tf_strat.pdf"),
  plot = p_long_med,
  width = 12,
  height = 4.5,
  units = "cm",
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_sd_med_ichorcna_tf_strat.pdf"))


# for each healthy sample, calcualte the total fraction between 100 and 155 bp
all_sd_healthy_100_155 <- all_sd_healthy %>%
  filter(x >= 100 & x <= 155) %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarize(total_fraction = sum(y)) %>%
  ungroup()

# for each cancer sample, calcualte the total fraction between 100 and 155 bp
all_sd_cancer_100_155 <- all_sd_cancer %>%
  filter(x >= 100 & x <= 155) %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarize(total_fraction = sum(y)) %>%
  ungroup()

# combine the two dataframes
all_sd_100_155 <- rbind(all_sd_healthy_100_155, all_sd_cancer_100_155)

# set the order of ichorcna_tf_strat, "Healthy" first
all_sd_100_155$ichorcna_tf_strat <- factor(all_sd_100_155$ichorcna_tf_strat,
  levels = ichorcna_tf_strat_hyper
)


# plot as box plot
p_box <- ggplot(all_sd_100_155, aes(x = ichorcna_tf_strat, y = total_fraction)) +
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
  # boxplot without outliers
  ylab("Sum of SD (100~155bp)") +
  xlab("ichorCNA Tumor Fraction") +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
  theme_classic() +
  # hide legend
  theme(legend.position = "none") +
  # increase the legend text size and line size
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.height = unit(0.9, "cm")
  ) +
  # scale_fill_nord(palette = "aurora") +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  # rotate x-axis tick and text by 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(
  filename = file.path(plot_saving_dir, "all_sd_100_155_boxplot.pdf"),
  plot = p_box,
  width = 15,
  height = 9,
  units = "cm",
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_sd_100_155_boxplot.pdf"))

# use patchwork to combine the  plots, hide all legends

p_long_med_nolegend <- p_long_med + theme(legend.position = "none")
p_box_nolegend <- p_box + theme(legend.position = "none")

ctdna_is_shorter <- p_long +
  (p_long_med_nolegend / p_box_nolegend) +
  plot_layout(widths = c(1.3, 1))

ggsave(
  filename = file.path(plot_saving_dir, "ctdna_len_is_more_variable_wide.pdf"),
  plot = ctdna_is_shorter,
  width = 18,
  height = 6,
  units = "cm",
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "ctdna_len_is_more_variable_wide.pdf"))
ggsave(
  filename = file.path(plot_saving_dir, "ctdna_len_is_more_variable_narrow.pdf"),
  plot = ctdna_is_shorter,
  width = 9,
  height = 5,
  units = "cm",
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "ctdna_len_is_more_variable_narrow.pdf"))





################################################################################
# plot p_long independently, width = 119.5 mm, height = 29.3 mm, save to pdf
################################################################################
p_long_new <- p_long +
  # axis text size
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5)) +
  # axis title size
  theme(axis.title.x = element_text(size = 5.5), axis.title.y = element_text(size = 5.5)) +
  # strip text size to 5
  theme(strip.text = element_text(size = 5))

ggsave(
  filename = file.path(plot_saving_dir, "all_sd_long_independent.pdf"),
  plot = p_long_new,
  width = 119.5 / 10,
  height = 69.1 / 20,
  units = "cm",
  dpi = 300
)

# plot p_long_med_nolegend independently, height = 60.4 mm, width = 81.4 mm, save to pdf

p_new <- (p_long_med_nolegend + p_box_nolegend) &

  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5)) &
  # axis title size
  theme(axis.title.x = element_text(size = 5.5), axis.title.y = element_text(size = 5.5)) &
  # remove x axis title
  theme(axis.title.x = element_blank())

ggsave(
  filename = file.path(plot_saving_dir, "all_sd_med_boxplot_independent.pdf"),
  plot = p_new,
  width =  90 / 10,
  height = 45 / 10,
  units = "cm",
  dpi = 300
)
