
library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)

plot_saving_dir <- "/home/nrlab/wang04/ulyses/plot_mae/plots/sd"

source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")

################################################################################
# sd
################################################################################
all_sd <- all_long_filter %>%
  dplyr::filter(assay == "sd") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)


# plot all samples together, color by author
((
  p_long <- ggplot(all_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.2, linewidth = 0.2) +
    ylab("SD") +
    xlab("Fragment Length (bp)") +
    theme_classic()
  
))

ggsave(
  filename = file.path(plot_saving_dir, "all_sd.pdf"),
  plot = p_long,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_sd.pdf"))


# plot median sd for each bicohor

all_sd_med <- all_sd %>%
  group_by(bicohort, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

((
  p_long_med <- ggplot(all_sd_med, aes(x, y, group = bicohort)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.4) +
    ylab("Median SD") +
    xlab("Fragment Length (bp)") +
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

# plot sd, facet by author
((
  p_long <- ggplot(all_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("SD") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # remove legend title
    guides(color = guide_legend(title = NULL)) +
    # set color palette
    scale_color_nord(palette = "aurora") +
    facet_wrap(~author)
))

ggsave(
  filename = file.path(plot_saving_dir, "all_sd_facet.pdf"),
  plot = p_long,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_sd_facet.pdf"))



# plot sd, facet by author and bicohort
((
  p_long <- ggplot(all_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("SD") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # hide legend
    #theme(legend.position = "none") +
    # remove legend title
    guides(color = guide_legend(title = NULL)) +
    # set color palette
    scale_color_manual(values = c("Cancer" = "grey", 
                                "Healthy" = "royalblue")) +
    facet_grid( bicohort ~ author)
))

ggsave(
  filename = file.path(plot_saving_dir, "all_sd_facet_author_cohort.pdf"),
  plot = p_long,
  width = 10,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_sd_facet_author_cohort.pdf"))


################################################################################
# plot sd facet by ichorcna_tf_strat
################################################################################
all_sd_healthy <- all_sd %>%
  filter(bicohort == "Healthy") %>%
  filter(ichorcna_tf_strat != 'unknown') %>%
  mutate(ichorcna_tf_strat = "Healthy")




all_sd_cancer <- all_sd %>%
  filter(bicohort == "Cancer") %>%
  filter(ichorcna_tf_strat != 'unknown')

# combine the two dataframes
all_sd_clean <- rbind(all_sd_healthy, all_sd_cancer)


label_levels <- c("Healthy", 
                  "(0, 0.01]", 
                  "(0.01, 0.03]", 
                  "(0.03, 0.1]", 
                  "(0.1, 0.2]", 
                  "(0.2, 1]")

# set the order of ichorcna_tf_strat, "Healthy" first
all_sd_clean$ichorcna_tf_strat <- factor(all_sd_clean$ichorcna_tf_strat, 
                                        levels = label_levels)

# count the number of sample for each ichorcna_tf_strat
make_label <- all_sd_clean %>%
  group_by(ichorcna_tf_strat) %>%
  summarize(n = length(unique(primary))) %>%
  ungroup() %>%
  mutate(facet_label = paste0("ichorCNA TF ", ichorcna_tf_strat, "\n", "(n=", n, ")"))

facet_label <- make_label$facet_label
names(facet_label) <- make_label$ichorcna_tf_strat


((
  p_long <- ggplot(all_sd_clean, aes(x, y, group = primary)) +
    #geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = Inf, fill = "lightgrey", alpha = 0.2)+
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("SD") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = c(50, 100, 167, 200, 250)) +
    # set color palette
    scale_color_nord(palette = "aurora") +
    facet_wrap(~ichorcna_tf_strat, labeller = labeller(ichorcna_tf_strat = facet_label)) +
    # set facet label backgroud color to light grey
    theme(strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"))

))


# save plot
ggsave(
  filename = file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.pdf"),
  plot = p_long,
  width = 7,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.pdf"))


################################################################################
# also use plotly to save p_long as html file
################################################################################
p_html <- ggplotly(p_long)

# save to htm
htmlwidgets::saveWidget(p_html, file.path(plot_saving_dir, "all_sd_facet_ichorcna_tf_strat.html"), selfcontained = TRUE)


################################################################################
# barplot: sum of sd for each ichorcna_tf_strat 
################################################################################
min <- 100
max <- 170

sum_sd_func <- function(input, min, max) {
  ans <- input %>%
    filter(x >= !!min & x <= !!max) %>%
    group_by(primary, ichorcna_tf_strat) %>%
    summarize(total_fraction = sum(y)) %>%
    ungroup()

  return(ans)
}

all_sd_healthy_sum_range <- sum_sd_func(all_sd_healthy, min, max)
all_sd_cancer_sum_range <- sum_sd_func(all_sd_cancer, min, max)

# combine the two dataframes
all_sd_sum_range <- rbind(all_sd_healthy_sum_range, all_sd_cancer_sum_range)

# set the order of ichorcna_tf_strat, "Healthy" first
all_sd_sum_range$ichorcna_tf_strat <- factor(all_sd_sum_range$ichorcna_tf_strat, 
                                        levels = label_levels)

# plot as box plot
  ylab_txt <- paste("Sum of SD (", min, "~", max, "bp fragments)", sep = "")
((
  p_box_sum_sd <- ggplot(all_sd_sum_range, aes(x = ichorcna_tf_strat, y = total_fraction)) +
  #boxplot without outliers
    geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA ) +
    geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.55) +
    ylab(ylab_txt) +
    xlab("ichorCNA Tumor Fraction") +
    # add p-value to the plot
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy") +
    theme_classic() +
    # remove legend
    theme(legend.position = "none") +
    scale_fill_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
))


ggsave(
  filename = file.path(plot_saving_dir, paste("all_sd_", min, "_", max, "_boxplot.pdf", sep = "")),
  plot = p_box_sum_sd,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, paste("all_sd_", min, "_", max, "_boxplot.pdf", sep = "")))

################################################################################
# barplot: sd of sd for each ichorcna_tf_strat
################################################################################

sd_sd_func <- function(input, min, max) {
  ans <- input %>%
    filter(x >= !!min & x <= !!max) %>%
    group_by(primary, ichorcna_tf_strat) %>%
    summarize(y_sd = sd(y)) %>%
    ungroup()
  return(ans)
}

all_sd_healthy_sd_range <- sd_sd_func(all_sd_healthy, min, max)
all_sd_cancer_sd_range <- sd_sd_func(all_sd_cancer, min, max)

# combine the two dataframes
all_sd_sd_range <- rbind(all_sd_healthy_sd_range, all_sd_cancer_sd_range)

# set the order of ichorcna_tf_strat, "Healthy" first
all_sd_sd_range$ichorcna_tf_strat <- factor(all_sd_sd_range$ichorcna_tf_strat, 
                                        levels = label_levels)

# plot as box plot
  sd_sd_ylab_txt <- paste("SD of SD (", min, "~", max, "bp fragments)", sep = "")
((
  p_box_sd_sd <- ggplot(all_sd_sd_range, aes(x = ichorcna_tf_strat, y = y_sd)) +
  #boxplot without outliers
    geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA ) +
    geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.55) +
    ylab(sd_sd_ylab_txt) +
    xlab("ichorCNA Tumor Fraction") +
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy") +
    theme_classic() +
    # remove legend
    theme(legend.position = "none") +
    scale_fill_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
))


ggsave(
  filename = file.path(plot_saving_dir, paste("all_sd_", min, "_", max, "_sd_sd_boxplot.pdf", sep = "")),
  plot = p_box_sd_sd,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, paste("all_sd_", min, "_", max, "_sd_sd_boxplot.pdf", sep = "")))

################################################################################
# median sd for each ichorcna_tf_strat
################################################################################

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
                                        levels = label_levels)

# plot median sd for each ichorcna_tf_strat
((
  p_long_med <- ggplot(all_sd_med, aes(x, y, group = ichorcna_tf_strat)) +
    geom_rect(xmin = min, xmax = max, ymin = 0, ymax = Inf, fill = "lightgrey", alpha = 0.2)+
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.95, linewidth = 0.7) +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    ylab("Median SD") +
    xlab("Fragment Length (bp)") +
    # add 155 and 167 bp to the x-axis tick and text
    scale_x_continuous(breaks = c(50, 80, 100, 155, 167, 200, 250)) +
    # only show x between 80 and 250
    coord_cartesian(xlim = c(80, 250)) +
    # remove legend title
    guides(color = guide_legend(title = NULL)) +
    # increase the legend text size and line size
    theme(legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8),
          legend.key.height = unit(0.9, "cm")) +
    theme_classic() +
        scale_color_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
))



ggsave(
  filename = file.path(plot_saving_dir, "all_sd_med_ichorcna_tf_strat.pdf"),
  plot = p_long_med,
  width = 6,
  height = 3,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_sd_med_ichorcna_tf_strat.pdf"))


################################################################################
# use patchwork to combine the  plots, hide all legends
################################################################################

p_long_med_nolegend <- p_long_med + theme(legend.position = "none")

ctdna_is_shorter <- p_long  + 
                    (p_long_med_nolegend / p_box_sum_sd) +
                    plot_layout(widths = c(1.3, 1))

ggsave(
  filename = file.path(plot_saving_dir, "ctdna_has_higher_sd.pdf"),
  plot = ctdna_is_shorter,
  width = 10,
  height = 6,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "ctdna_has_higher_sd.pdf"))





