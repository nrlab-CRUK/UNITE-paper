# libs and sources
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(MultiAssayExperiment))
suppressPackageStartupMessages(library(nord))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(kableExtra))


source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")
source("/home/nrlab/wang04/ulyses/plot_mae/viz_hyperparmeters.R")

# parameters
# mae <- readRDS(mae_file_latest)

# make plot_saving_dir with iteration_num_hyper
plot_saving_dir_base <- "/home/nrlab/wang04/ulyses/plot_mae/plots/sample_alluvial"
plot_saving_dir <- file.path(plot_saving_dir_base, paste0("iteration_", iteration_num_hyper))
if (!dir.exists(plot_saving_dir)) {
  dir.create(plot_saving_dir)
}


# meta <- colData(mae)

meta <- read_csv(meta_csv_latest) %>%
  # only keep plasma
  filter(sample_type %in% sample_type_hyper) %>%
  # remove the Stephen Cristiano samples
  filter(author != "Stephen Cristiano") %>%
  # keep authors in authors_hyper
  filter(author %in% authors_hyper) %>%
  # rename bam_id to primary
  rename(primary = bam_id)



# meta QC and tracking problems

meta %>%
  filter(is.na(seq_depth)) |>
  write_csv("/home/nrlab/wang04/ulyses/plot_mae/plots/sample_alluvial/samples_without_seq_depth.csv")

# get the column data
mae_raw <- meta

# figure out cohort less than 10
cohort_less_than_20 <- mae_raw %>%
  group_by(cohort) %>%
  summarise(n = n()) %>%
  filter(n < 20) %>%
  pull(cohort)

# change cohort to "other" if cohort is less than 20
mae_raw <- mae_raw %>%
  mutate(cohort2 = case_when(
    cohort %in% cohort_less_than_20 ~ "other",
    TRUE ~ cohort
  ))

# add `outlier` column to mae raw
message("Adding outlier")
mae_raw <- mae_raw %>%
  mutate(outlier = case_when(
    primary %in% outliers_hyper ~ TRUE,
    TRUE ~ FALSE
  ))

# add tri-cohort column to mae raw
message("Adding tricohort")
mae_raw <- mae_raw %>%
  mutate(tricohort = case_when(
    cohort == "Healthy" ~ "Healthy",
    cohort %in% noncancer_disease_hyper ~ "NCD",
    TRUE ~ "Cancer"
  ))
# add tri-timepoint column to mae raw
message("Adding tri-timepoint")

mae_raw <- mae_raw %>%
  mutate(tritimepoint = case_when(
    timepoint == "1" ~ "Earliest",
    is.na(timepoint) ~ "Unknown",
    TRUE ~ "Later"
  ))

# add ichorcna_tf_strat to mae raw
message("Adding ichorcna_tf_strat")
mae_raw_h <- mae_raw %>%
  filter(tricohort == "Healthy") %>%
  mutate(ichorcna_tf_strat = "Healthy")

# mae_raw_ncd <- mae_raw %>%
#  filter(tricohort == "NCD") %>%
#  mutate(ichorcna_tf_strat = "NCD")

mae_raw_c <- mae_raw %>%
  filter(tricohort %in% c("Cancer", "NCD"))

mae_raw_c <- mae_raw_c %>%
  add_ichorcna_tf_strat_col()

# mutate(ichorcna_tf_strat = case_when(
#   ichorcna_tf <= 0.03 ~ "[0, 0.03]",
#   # ichorcna_tf <= 0.03 & ichorcna_tf > 0.01 ~ "(0.01, 0.03]",
#   ichorcna_tf <= 0.1 & ichorcna_tf > 0.03 ~ "(0.03, 0.1]",
#   ichorcna_tf <= 0.2 & ichorcna_tf > 0.1 ~ "(0.1, 0.2]",
#   ichorcna_tf <= 1 & ichorcna_tf > 0.2 ~ "(0.2, 1]",
#   TRUE ~ "unknown"
# ))

mae_raw_tidy <- rbind(
  mae_raw_h,
  # mae_raw_ncd,
  mae_raw_c
)

# set author factor level to authos_hyper



# add a column 'selected' to mae raw
mae_raw_tidy <- mae_raw_tidy %>%
  mutate(selected = case_when(
    primary %in% mae_col_data$primary ~ "Yes",
    TRUE ~ "No"
  ))


# summary of the count of each author, tricohort
mae_raw_tidy %>%
  group_by(author, tricohort) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  arrange(author) %>%
  write_csv("/home/nrlab/wang04/ulyses/plot_mae/plots/sample_alluvial/QC/author_tricohort.csv")

# plot as alluvial plot using ggplot, color by `selected`
strata <- c(
  "seq_depth",
  "selected",
  "data_source",
  "author",
  "tricohort",
  "cohort2",
  "timepoint",
  "tritimepoint",
  "ichorcna_tf_strat",
  "outlier",
  "stage",
  "library_kit",
  "extraction_kit",
  "seq_platform",
  "sample_type"
)

alluvial_data <- select(mae_raw_tidy, all_of(strata))

# stratefy the seq_depth column into 2 groups: <0.1x and >=0.1x
alluvial_data <- alluvial_data %>%
  mutate(tri_seq_depth = case_when(
    seq_depth < 0.1 ~ "<0.1x",
    is.na(seq_depth) ~ "NA",
    TRUE ~ ">=0.1x"
  ))


# make the NA in all columns as "NA"
alluvial_data <- alluvial_data %>%
  mutate_all(~ ifelse(is.na(.), "NA", .))

# set the factor levels of data_source
alluvial_data$data_source <- factor(alluvial_data$data_source,
  levels = data_source_hyper
)

# set the factor levels of timepoint

alluvial_data$timepoint <- factor(alluvial_data$timepoint,
  levels = tp_levels_hyper
)
# do the alluvial plot

strata_changed <- c(
  "tri_seq_depth",
  "selected",
  "data_source",
  "author",
  "tricohort",
  "cohort2",
  # "timepoint",
  "tritimepoint",
  "ichorcna_tf_strat",
  "outlier"
)

# rename the author name
alluvial_data$author <- author_rename_vec_hyper[alluvial_data$author]


alluvial_data_wide <- alluvial_data %>%
  group_by_at(strata_changed) %>%
  summarise(freq = n()) %>%
  ungroup()

# set the factor levels of ichorcna_tf_strat
alluvial_data_wide$ichorcna_tf_strat <- factor(alluvial_data_wide$ichorcna_tf_strat,
  levels = ichorcna_tf_strat_col_levels_alluvial_hyper
)



# remoe the 'et al' from author
# alluvial_data_wide$author <- gsub(" et al", "", alluvial_data_wide$author)



base::options(ggalluvial.decreasing = FALSE)

n_sample_total <- sum(alluvial_data_wide$freq)
n_sample_selected <- sum(alluvial_data_wide$freq[alluvial_data_wide$selected == "Yes"])

############################################################################################################
# alluvial
############################################################################################################

# save alluvial_data_wide to csv file
alluvial_data_wide %>%
  write_csv(file.path(plot_saving_dir, "alluvial_wide_data.csv"))

message("saved to", file.path(plot_saving_dir, "alluvial_wide_data.csv"))

n_colors <- alluvial_data_wide$author |>
  unique() |>
  length()

library(viridis)
# color_values <- viridis(length(authors))
# color_values <- inferno(length(authors))
color_values <- plasma(length(authors_hyper))
# color_values <- colorspace::diverge_hcl(length(authors))
# color_values <- turbo(length(authors))



p2 <- ggplot(
  alluvial_data_wide,
  aes(
    y = freq,
    # axis1 = author,
    # axis2 = data_source,
    # axis3 = tri_seq_depth,
    # axis4 = tricohort,
    # axis5 = ichorcna_tf_strat,
    # axis6 = tritimepoint,
    # axis7 = selected
    axis1 = author,
    axis2 = data_source,
    # axis3 = ichorcna_tf_strat,
    axis3 = tri_seq_depth,
    axis4 = tricohort,
    axis5 = tritimepoint,
    axis6 = selected
  )
) +
  geom_alluvium(aes(fill = `author`), width = 1 / 7, alpha = 0.98) +
  # geom_flow(aes(fill = author), width = 0.2, alpha = 1) +
  geom_stratum(
    fill = "black",
    size = 0.1,
    width = 1 / 7,
    color = "#ffffff",
    alpha = 0.61
  ) +
  # geom_lode() +
  # geom_label(stat = "stratum", aes(label = after_stat(stratum)), size = 0.4, label.size = 0, angle = -90) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "white", angle = 0, size = 2) +
  # ggrepel::geom_text_repel(aes(label = author), stat = "stratum", size = 4, direction = "y", nudge_x = -.5) +
  scale_x_discrete(expand = c(0, 0), limits = c("Author", "Database", "Coverage", "Cohort", "Timepoints", "Selected")) +
  # scale_fill_manual(values = rev(nord_palettes$afternoon_prarie)) +
  # use color pelate from RColorBrewer
  # scale_fill_brewer(type = "div") +
  # use colorspace::diverge_hcl() to make a diverging color palette
  scale_fill_manual(values = colorspace::diverge_hcl(n_colors)) +
  # add x axis stratums labels
  labs(x = NULL, y = NULL) +
  # add y axis to both left and right
  scale_y_continuous(expand = c(0, 0), position = "right", breaks = c(0, n_sample_selected, n_sample_total)) +
  # add geom_hline at total number of samples
  # geom_hline(yintercept = n_sample_total, color = "black", linewidth = 0.2, linetype = "dashed") +
  # add geom_hline at total number of selected samples (yes)
  # geom_hline(yintercept = n_sample_selected, color = "black", linewidth = 0.2, linetype = "dashed") +
  # y axis to the right side of the plot
  # coord_flip() +
  theme_classic() +
  theme(
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    # axis axis text size to 5
    axis.text = element_text(size = 5.5)
  ) +
  # make y axis thicker
  theme(axis.line.y = element_line(linewidth = 0.8)) +
  theme(legend.position = "none") +
  # remove any plot margin
  theme(plot.margin = unit(c(0, 0, 0, 10), "mm"))

ggsave(filename = file.path(plot_saving_dir, "alluvial.pdf"), p2, width = 165, height = 90, units = "mm", dpi = 300)
message(file.path(plot_saving_dir, "alluvial.pdf"))


############################################################################################################
# plot pie chart of author, data_source, tricohort, cohort, timepoint, ichorcna_tf_strat, selected
############################################################################################################


# a function for plotting the stacked bar plot of each variables in alluvial plot

alluvial_bar_plot <- function(variable,
                              alluvial_data,
                              plot_saving_dir = getwd(),
                              plot_width = 45,
                              plot_height = 52,
                              plot_unit = "mm",
                              plot_dpi = 300) {
  # plot pie chart of author
  variable_data <- alluvial_data %>%
    select({{ variable }}, selected) %>%
    group_by(pick(everything())) %>%
    summarise(freq = n())

  print(variable_data, Inf)

  # set the levels of selected
  variable_data$selected <- factor(variable_data$selected, levels = c("Yes", "No"))

  # stack bar plot, x axis is variable, y is frequency, color by selected, stack by selected
  p <- ggplot(variable_data, aes(x = .data[[variable]], y = freq, fill = selected)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = variable, y = "Frequency") +
    # add labels
    theme_classic() +
    # remove legend
    theme(legend.position = "none") +
    # add labells to the selected==yes bar
    geom_text(aes(label = freq), position = position_stack(vjust = 0.5), color = "black", size = 1.8) +
    # rotate x axis text 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # manual color
    scale_fill_manual(values = c(selected_yes_color, selected_no_color)) +
    theme(legend.position = "none") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
    theme(axis.title = element_text(size = 5.5)) +
    theme(axis.text = element_text(size = 5))

  plot_file_name <- paste0(variable, "_bar.pdf")
  plotfilename <- file.path(plot_saving_dir, plot_file_name)

  # get the directory of the filename
  dirname(plotfilename) %>%
    dir.create(recursive = TRUE, showWarnings = FALSE)

  ggsave(
    filename = plotfilename,
    p,
    width = plot_width,
    height = plot_height,
    units = plot_unit,
    dpi = plot_dpi
  )

  message("Saved plot to ", plotfilename)
  return(p)
}




variable_list <- c(
  "author",
  # "seq_depth",
  "tri_seq_depth",
  "data_source",
  "tricohort",
  "cohort2",
  "timepoint",
  "ichorcna_tf_strat",
  "selected"
)

tmp <- lapply(variable_list,
  alluvial_bar_plot,
  alluvial_data = alluvial_data,
  plot_saving_dir = file.path(plot_saving_dir, "QC"),
  plot_width = 70,
  plot_height = 60,
  plot_unit = "mm",
  plot_dpi = 300
)

names(tmp) <- variable_list


# plot seq depth
seq_depth_data <- alluvial_data %>%
  select(seq_depth, selected) %>%
  mutate(seq_depth = as.numeric(seq_depth))

# set the levels of selected
seq_depth_data$selected <- factor(seq_depth_data$selected, levels = c("Yes", "No"))

# plot as stacked bar plot, color by selected
p_depth <- ggplot(seq_depth_data, aes(x = seq_depth, fill = selected)) +
  # histogram
  geom_histogram(binwidth = 1) +
  labs(x = "Sequencing Depth", y = "Frequency", fill = "Dataset") +
  # add text labels to each stacked bar regarding count number
  # stat_bin(geom = "text", aes(label = ..count..), vjust = 0.5, size = 1, position = "stack") +
  theme_classic() +
  # make x axis text size smaller
  # color by aurora color palete from nord package
  scale_fill_manual(values = c(selected_yes_color, selected_no_color)) +
  theme(legend.position = "none") +
  # rotate x axis 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # axis text size = 5.5
  theme(axis.text = element_text(size = 5)) +
  # axis title size = 5.5
  theme(axis.title = element_text(size = 5.5))

# ggbreak break y axis between 500 and 1400
p_depth2 <- p_depth +
  ggbreak::scale_y_break(c(500, 1500), ticklabels = c(1500, 1570), expand = TRUE)

p_depth2 <- print(p_depth2)
# zoom
# p_zoom <- p_depth + coord_cartesian(xlim = c(0, 10))

ggsave(filename = file.path(plot_saving_dir, "seq_depth_histo.pdf"), p_depth, width = 76, height = 40, units = "mm", dpi = 300)
ggsave(filename = file.path(plot_saving_dir, "seq_depth_histo_break_y.pdf"), p_depth2, width =76, height = 40, units = "mm", dpi = 300)
# ggsave(filename = file.path(plot_saving_dir, "seq_depth_histo_zoom.pdf"), p_zoom, width = 70, height = 60, units = "mm", dpi = 300)

message(file.path(plot_saving_dir, "seq_depth_histo.pdf"))
message(file.path(plot_saving_dir, "seq_depth_histo_break_y.pdf"))
# message(file.path(plot_saving_dir, "seq_depth_histo_zoom.pdf"))


# arrange the qc plots using patchwork
((

  patchwork1 <- tmp$author + tmp$cohort2 + tmp$data_source + tmp$tricohort +
    plot_layout(ncol = 4, nrow = 1, widths = c(1, 1, 0.3, 0.3))
))
((

patchwork2 <- tmp$timepoint + p_depth + tmp$selected +
  plot_layout(ncol = 3, nrow = 1, widths = c(1, 1.3, 0.3))

))

final <- patchwork1 / patchwork2
ggsave(
  filename = file.path(plot_saving_dir, "alluvial_qc.pdf"),
  final,
  width = 180,
  height = 105,
  units = "mm",
  dpi = 300
)
message(file.path(plot_saving_dir, "alluvial_qc.pdf"))




############################################################################################################
# do a stacked bar plot colored by author, x axis is cohort, y is frequency
############################################################################################################


alluvial_data_selected <- mae_raw_tidy %>%
  filter(selected == "Yes" | selected == TRUE)

cohort_levels <- group_by(alluvial_data_selected, cohort) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  pull(cohort)

# set the factor levels of cohort
alluvial_data_selected$cohort <- factor(alluvial_data_selected$cohort, levels = cohort_levels)
# sort the cohort by frequency

# set the factor levels of author
alluvial_data_selected$author <- factor(alluvial_data_selected$author, levels = authors_hyper)


p3 <- ggplot(alluvial_data_selected, aes(x = cohort, fill = author)) +
  geom_bar() +
  labs(x = "Cohort", y = "Frequency") +
  theme_classic() +
  # color by aurora color palete from nord package
  scale_fill_manual(values = author_color_values_hyper) +
  # rotate x axis 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")

ggsave(
  filename = file.path(plot_saving_dir, "stacked_bar.pdf"),
  plot = p3,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)

message("file saved to ", file.path(plot_saving_dir, "stacked_bar.pdf"))
############################################################################################################
# do a stacked bar plot colored by author, x axis is cohort, y is frequency
############################################################################################################

alluvial_data_selected <- mae_raw_tidy %>%
  filter(selected == TRUE | selected == "Yes")

cohort_levels <- group_by(alluvial_data_selected, cohort) %>%
  summarise(n = n()) %>%
  # sort the cohort by frequency
  arrange(desc(n)) %>%
  pull(cohort)

alluvial_data_selected$cohort <- factor(alluvial_data_selected$cohort, levels = cohort_levels)
alluvial_data_selected$author <- factor(alluvial_data_selected$author, levels = authors_hyper)

p3 <- ggplot(alluvial_data_selected, aes(x = cohort, fill = author)) +
  geom_bar() +
  labs(x = "Cohort", y = "Frequency", fill = "Dataset") +
  theme_classic() +
  # make x axis text size smaller
  theme(axis.text.x = element_text(size = 7)) +
  # color by aurora color palete from nord package
  scale_fill_manual(values = author_color_values_hyper_viridis) +
  # rotate x axis 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")

ggsave(
  filename = file.path(plot_saving_dir, "all_selected_stacked_bar.pdf"), p3,
  width = 180,
  height = 90,
  units = "mm",
  dpi = 300
)
message(file.path(plot_saving_dir, "all_selected_stacked_bar.pdf"))
############################################################################################################
# plot for each ichorcna_tf_strat
############################################################################################################
# write a function for plotting stacked bar plot

stack_bar_f <- function(x, dir, height = 45, width = 90, units = "mm", dpi = 300, legend = TRUE) {
  if (legend == FALSE) {
    plot_file_name <- paste0("stacked_bar_", unique(x$ichorcna_tf_strat), "_no_legend.pdf")
  } else {
    plot_file_name <- paste0("stacked_bar_", unique(x$ichorcna_tf_strat), ".pdf")
  }
  plot_file_name_full <- file.path(dir, plot_file_name)

  cohort_levels <- group_by(x, cohort) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    pull(cohort)

  # set the factor levels of cohort
  x$cohort <- factor(x$cohort, levels = cohort_levels)

  max_y <- x %>%
    group_by(cohort) %>%
    summarize(n = n()) %>%
    pull("n") %>%
    max()
  y_label <- unique(x$ichorcna_tf_strat)
  y_label <- paste0(y_label, " (n)")

  p <- ggplot(x, aes(x = cohort)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.2, size = 2) +
    geom_bar(aes(fill = author)) +
    # expand_limits, set y from 0 to max(y) + 0.02 * max(y)
    expand_limits(y = c(0, max_y + 0.1 * max_y)) +
    labs(
      x = unique(x$ichorcna_tf_strat),
      y = y_label,
      fill = "Dataset"
      # subtitle = unique(x$ichorcna_tf_strat)
    ) +
    # remove x axis
    theme_classic() +
    # make x axis text size smaller
    theme(axis.text = element_text(size = 5)) +
    # axis title size = 6
    theme(axis.title = element_text(size = 5.5)) +
    # make legend size smaller
    theme(legend.text = element_text(size = 5.5)) +
    # subtitle size to 5.5
    theme(plot.subtitle = element_text(size = 6)) +
    # color by aurora color palete from nord package
    scale_fill_manual(values = author_color_values_hyper) +
    # rotate x axis 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "right") +
    theme(axis.title.x = element_blank())

  if (legend == FALSE) {
    p <- p + theme(legend.position = "none")
  }

  ggsave(
    filename = plot_file_name_full,
    p,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )

  message("Saved plot to ", plot_file_name_full)

  return(p)
}

alluvial_data_selected <- mae_raw_tidy %>%
  filter(selected == TRUE | selected == "Yes") %>%
  filter(cohort != "Healthy")

# set author factor level to authors_hyper
alluvial_data_selected$author <- factor(alluvial_data_selected$author, levels = authors_hyper)

alluvial_data_selected_list <- alluvial_data_selected %>%
  group_by(ichorcna_tf_strat) %>%
  group_split()

tmp <- lapply(alluvial_data_selected_list,
  stack_bar_f,
  dir = plot_saving_dir,
  height = 70,
  width = 90,
  units = "mm",
  dpi = 300
)

tmp <- lapply(alluvial_data_selected_list,
  stack_bar_f,
  dir = plot_saving_dir,
  height = 37,
  width = 85.3,
  units = "mm",
  dpi = 300,
  legend = FALSE
)





################################################################################
# plot the stacked bar plot for each ichorcna_tf_strat
################################################################################



plot_file_name <- paste0("stacked_bar_", "ichorcna_TF_strat", ".pdf")
plot_file_name_full <- file.path(plot_saving_dir, plot_file_name)
plot_file_name_no_legend <- paste0("stacked_bar_", "ichorcna_TF_strat_no_legend", ".pdf")
plot_file_name_full_no_legend <- file.path(plot_saving_dir, plot_file_name_no_legend)

alluvial_data_selected <- mae_raw_tidy %>%
  filter(selected == TRUE | selected == "Yes")


alluvial_data_selected <- alluvial_data_selected %>%
  # set "Healthy" as the first level of ichorcna_tf_strat
  mutate(ichorcna_tf_strat = factor(ichorcna_tf_strat, levels = ichorcna_tf_strat_col_levels_hyper))

# set author factor level to authors_hyper
alluvial_data_selected$author <- factor(alluvial_data_selected$author, levels = authors_hyper)

# with legend
p <- ggplot(alluvial_data_selected, aes(x = ichorcna_tf_strat)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.2, size = 2) +
  geom_bar(aes(fill = author)) +
  labs(x = "ichorCNA TF Strata", y = "Count (n)", fill = "Dataset") +
  # add count number to each bar
  theme_classic() +
  # make x axis text size smaller
  theme(axis.text = element_text(size = 5.5)) +
  # make legend size smaller
  theme(legend.text = element_text(size = 5)) +
  # color by aurora color palete from nord package
  scale_fill_manual(values = author_color_values_hyper) +
  # rotate x axis 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "right")


ggsave(
  filename = plot_file_name_full,
  p,
  height = 60,
  width = 90,
  units = "mm",
  dpi = 300
)

message("Saved plot to ", plot_file_name_full)

# without legend
p_no_legend <- p +
  # make y between 0 and 800
  coord_cartesian(ylim = c(0, 800)) +
  theme(legend.position = "none") +
  # set y-axis title as "Count (n)"
  labs(y = "Count (n)") +
  # axis title font size = 6
  theme(axis.title = element_text(size = 5.5)) +
  # axis text font size = 6
  theme(axis.text = element_text(size = 5.5))

ggsave(
  filename = plot_file_name_full_no_legend,
  p_no_legend,
  height = 53,
  width = 34,
  units = "mm", dpi = 300
)

message("Saved plot to ", plot_file_name_full_no_legend)


################################################################################
# plot box plot for ichorcna_tf, x is ichorcna_tf_strat, y is ichorcna_tf
################################################################################

p_ichorcna_tf <- ggplot(alluvial_data_selected, aes(x = ichorcna_tf_strat, y = ichorcna_tf)) +
  geom_violin(fill = "lightgrey", alpha = 0.8, linewidth = 0.1, width = 1, color = "lightgrey") +
  geom_boxplot(aes(fill = ichorcna_tf_strat), width = 0.1, alpha = 0.7, linewidth = 0.1, outlier.shape = NA) +
  # add horizontal dashed lines at 0.03, 0.1
  geom_hline(yintercept = 0.03, color = "black", linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = 0.1, color = "black", linetype = "dashed", size = 0.1) +
  # add y axis breaks to 0, 0.03, 0.1, 0.2, 0.4, 0.6
  scale_y_continuous(breaks = c(0, 0.03, 0.1, 0.2, 0.4, 0.6)) +
  # geom_jitter(aes(fill = ichorcna_tf_strat),
  #   color = "black",
  #   width = 0.25,
  #   alpha = 0.5,
  #   size = 0.25,
  #   shape = 21,
  #   stroke=0.1) +
  # geom_boxplot(aes(color = "darkgrey")) +
  labs(x = "ichorCNA TF Strata", y = "Tumor Fraction inferred by ichorCNA") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # axis title text size = 5.5
  theme(axis.title = element_text(size = 5.5)) +
  # axis text size = 5.5
  theme(axis.text = element_text(size = 5.5)) +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  theme(legend.position = "none") +
  # remove x axis title
  theme(axis.title.y = element_blank()) +
  # flip the x and y axis
  coord_flip()

ggsave(
  filename = file.path(plot_saving_dir, "ichorcna_tf_boxplot.pdf"),
  p_ichorcna_tf,
  height = 53,
  width = 66,
  units = "mm",
  dpi = 300
)
message(file.path(plot_saving_dir, "ichorcna_tf_boxplot.pdf"))


# combine strat and tf plot using patchwork
p_combined <- p_ichorcna_tf +
  p_no_legend +
  plot_layout(widths = c(1, 0.32))

ggsave(
  filename = file.path(plot_saving_dir, "ichorcna_tf_combined.pdf"),
  p_combined, height = 53, width = 95, units = "mm", dpi = 300
)
message(file.path(plot_saving_dir, "ichorcna_tf_combined.pdf"))


# only plot the boxplot for ichorcna_tf_strat when ichorcna_tf_strat is in c("Healthy", "[0, 0.03]", "(0.03, 0.1]")
alluvial_data_selected_zoom <- alluvial_data_selected %>%
  filter(ichorcna_tf_strat %in% c("Healthy", "[0, 0.03]"))

p_ichorcna_tf_zoom <- ggplot(alluvial_data_selected_zoom, aes(x = ichorcna_tf_strat, y = ichorcna_tf)) +
  geom_violin(fill = "lightgrey", alpha = 0.8, linewidth = 0.1, width = 1, color = "lightgrey") +
  geom_boxplot(aes(fill = ichorcna_tf_strat), width = 0.1, alpha = 0.7, linewidth = 0.1, outlier.shape = NA) +
  # add horizontal dashed lines at 0.03, 0.1
  geom_hline(yintercept = 0.03, color = "black", linetype = "dashed", size = 0.1) +
  # geom_hline(yintercept = 0.1, color = "black", linetype = "dashed", size = 0.1) +
  # add y axis breaks to 0, 0.03, 0.1, 0.2, 0.4, 0.6
  scale_y_continuous(breaks = c(0, 0.03)) +
  geom_jitter(aes(fill = ichorcna_tf_strat),
    color = "black",
    width = 0.25,
    alpha = 0.5,
    size = 0.2,
    shape = 21,
    stroke = 0.02
  ) +
  labs(x = "ichorCNA TF Strata", y = "Tumor Fraction inferred by ichorCNA") +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # axis title text size = 5.5
  theme(axis.title = element_text(size = 5)) +
  # axis text size = 5.5
  theme(axis.text = element_text(size = 5)) +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  theme(legend.position = "none") +
  # remove x axis title
  theme(axis.title = element_blank()) +
  # flip the x and y axis
  coord_flip()

ggsave(
  filename = file.path(plot_saving_dir, "ichorcna_tf_plot_zoom_in.pdf"),
  p_ichorcna_tf_zoom,
  height = 25,
  width = 30,
  units = "mm",
  dpi = 300
)

message(file.path(plot_saving_dir, "ichorcna_tf_plot_zoom_in.pdf"))





################################################################################
# scatterpie plot of cohort and author
# https://cran.r-project.org/web/packages/scatterpie/vignettes/scatterpie.html
# https://stackoverflow.com/questions/56518292/how-to-plot-a-scatterpie-plot-with-ggplot2


# get the column data
mae_raw <- colData(mae) %>%
  as_tibble(rownames = "primary") %>%
  filter(author != "Stephen Cristiano")

# add `outlier` column to mae raw
message("Adding outlier")
mae_raw <- mae_raw %>%
  mutate(outlier = case_when(
    primary %in% outliers_hyper ~ TRUE,
    TRUE ~ FALSE
  ))

# add tri-cohort column to mae raw
message("Adding tricohort")
mae_raw <- mae_raw %>%
  mutate(tricohort = case_when(
    cohort == "Healthy" ~ "Healthy",
    cohort %in% noncancer_disease_hyper ~ "NCD",
    TRUE ~ "Cancer"
  ))

# add ichorcna_tf_strat to mae raw
message("Adding ichorcna_tf_strat")
mae_raw_h <- mae_raw %>%
  filter(tricohort == "Healthy") %>%
  mutate(ichorcna_tf_strat = "Healthy")

mae_raw_ncd <- mae_raw %>%
  filter(tricohort == "NCD") %>%
  mutate(ichorcna_tf_strat = "NCD")

mae_raw_c <- mae_raw %>%
  filter(tricohort == "Cancer")

# add ichorcna_tf_strat to mae_raw_c
mae_raw_c <- add_ichorcna_tf_strat_col(mae_raw_c)



# mutate(ichorcna_tf_strat = case_when(
#   ichorcna_tf <= 0.01 ~ "[0, 0.01]",
#   ichorcna_tf <= 0.03 & ichorcna_tf > 0.01 ~ "(0.01, 0.03]",
#   ichorcna_tf <= 0.1 & ichorcna_tf > 0.03 ~ "(0.03, 0.1]",
#   ichorcna_tf <= 0.2 & ichorcna_tf > 0.1 ~ "(0.1, 0.2]",
#   ichorcna_tf <= 1 & ichorcna_tf > 0.2 ~ "(0.2, 1]",
#   TRUE ~ "Unknown"
# ))



mae_raw_tidy <- rbind(mae_raw_h, mae_raw_ncd, mae_raw_c)




# add a column 'selected' to mae raw
mae_raw_tidy <- mae_raw_tidy %>%
  mutate(selected = case_when(
    primary %in% mae_col_data$primary ~ TRUE,
    TRUE ~ FALSE
  ))



# plot as alluvial plot using ggplot, color by `selected`
strata <- c(
  "selected",
  "data_source",
  "author",
  "tricohort",
  "cohort",
  "timepoint",
  "ichorcna_tf_strat",
  "outlier",
  "stage",
  "library_kit",
  "extraction_kit",
  "seq_platform",
  "sample_type"
)

alluvial_data <- select(mae_raw_tidy, all_of(strata))

# make the NA in all columns as "NA"
alluvial_data <- alluvial_data %>%
  mutate_all(~ ifelse(is.na(.), "NA", .))


# find the cohort with less than 10 samples
cohort_less_than_10 <- alluvial_data %>%
  filter(selected == TRUE) %>%
  group_by(cohort) %>%
  summarise(n = n()) %>%
  filter(n < 10) %>%
  pull(cohort)


tibble_other <- alluvial_data %>%
  filter(selected == TRUE | selected == "Yes") %>%
  group_by(cohort) %>%
  summarise(n = n()) %>%
  filter(n < 10)

# change cohort to "other" if cohort is less than 10
alluvial_data <- alluvial_data %>%
  mutate(cohort2 = case_when(
    cohort %in% cohort_less_than_10 ~ "other",
    TRUE ~ cohort
  ))


library(scatterpie)

scatterpie_data <- alluvial_data %>%
  filter(selected == TRUE | selected == "Yes") %>%
  select(cohort2, author)

cohorts <- scatterpie_data$cohort2 %>% unique()
authors <- scatterpie_data$author %>% unique()

sp_data2 <- scatterpie_data %>%
  group_by(cohort2, author) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  # pivot_wider to make the data wide, names from author, values from n
  pivot_wider(names_from = author, values_from = n)


# make a grid index of x and y, 5 rows, 6 columns

pos <- expand.grid(long = 1:5, lat = 1:6)[1:nrow(sp_data2), ]


sp_data3 <- sp_data2 %>%
  # make two vectors for the x and y position of the pie chart
  mutate(
    long = pos$long,
    lat = pos$lat,
    total = rowSums(sp_data2[, 2:ncol(sp_data2)], , na.rm = TRUE)
  ) %>%
  # change all NA to 0
  mutate_at(vars(-cohort2, -long, -lat), ~ replace_na(., 0))
# set the factor levels of author


# scatterpie plot

p <- ggplot() +
  geom_scatterpie(
    data = sp_data3,
    aes(x = long, y = lat, r = 0.0015 * (total)),
    alpha = 0.99,
    cols = sort(authors),
    color = NA
  ) +
  coord_equal() +
  # add number of samples and cohort name below each pie chart
  geom_text(
    data = sp_data3,
    aes(x = long, y = lat - 0.2, label = paste0("n=", total)),
    size = 3,
    color = "black"
  ) +
  # add cohort name
  geom_text(
    data = sp_data3,
    aes(x = long, y = lat + 0.2, label = cohort2),
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(values = author_color_values_hyper) +
  theme_void()

#



ggsave(
  filename = file.path(plot_saving_dir, "scatterpie.pdf"), p,
  width = 180, height = 140, units = "mm", dpi = 300
)

message("Saved plot to ", file.path(plot_saving_dir, "scatterpie.pdf"))




# table_obj <- kable(tibble_other, "html") %>%
#   kable_styling()

other_cohort_csv_file <- file.path(plot_saving_dir, "other_cohort_details.csv")
write_csv(tibble_other, other_cohort_csv_file)
message("Saved Other cohort table to ", other_cohort_csv_file)
