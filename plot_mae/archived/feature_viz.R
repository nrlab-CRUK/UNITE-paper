
library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)

mae <- readRDS("/scratchb/nrlab/wang04/ulyses_feat_viz/MAE3.rds")

plot_saving_dir <- "/home/nrlab/wang04/ulyses/plot_mae/plots/length"



# healthy -----------------------------------------------------------------

################################################################################
# length
################################################################################


all_h <- mae[, mae$cohort == 'Healthy', c("sd", "length")]

all_h_long <- longFormat(all_h, colDataCols = "author") %>%
  tibble::as_tibble() 

all_h_len<- all_h_long %>%
  dplyr::filter(assay == "length") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)


# plot all samples together, color by author
((
  p_long <- ggplot(all_h_len, aes(x, y, group = primary)) +
    geom_line(aes(color = author), alpha = 0.2, linewidth = 0.2) +
    ylab("Fraction") +
    xlab("Fragment Length (bp)") +
    theme_classic()
  
))

ggsave(
  filename = file.path(plot_saving_dir, "all_h_len.pdf"),
  plot = p_long,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_h_len.pdf"))

# calculate the median length fraction for each sample, color by author
all_h_len_med <- all_h_len %>%
  group_by(author, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

((
  p_long_med <- ggplot(all_h_len_med, aes(x, y, group = author)) +
    geom_line(aes(color = author), alpha = 0.8, linewidth = 0.4) +
    ylab("Fraction") +
    xlab("Fragment Length (bp)") +
    theme_classic()
  
))

ggsave(
  filename = file.path(plot_saving_dir, "all_h_len_med.pdf"),
  plot = p_long_med,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_h_len_med.pdf"))

# plot length, facet by author
((
  p_long <- ggplot(all_h_len, aes(x, y, group = primary)) +
    geom_line(aes(color = author), alpha = 0.3, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("Fraction") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # set color palette
    scale_color_nord(palette = "aurora") +
    facet_wrap(~author)
))

ggsave(
  filename = file.path(plot_saving_dir, "all_h_len_facet.pdf"),
  plot = p_long,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_h_len_facet.pdf"))

################################################################################
# SD
################################################################################

# plot sd, facet by author
all_h_sd<- all_h_long %>%
  dplyr::filter(assay == "sd") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)

((
  p_sd <- ggplot(all_h_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = author), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("SD of length across 5MB bins") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # set color palette
    scale_color_nord(palette = "aurora") +
    facet_wrap(~author)
))

ggsave(
  filename = file.path(plot_saving_dir, "all_h_sd_facet.pdf"),
  plot = p_sd,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_h_sd_facet.pdf"))



# plot everything together

((
  p_sd_all <- ggplot(all_h_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = author), alpha = 0.7, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("SD of length across 5MB bins") +
    theme_classic() +
    # hide legend
    # set color palette
    scale_color_nord(palette = "aurora")
))

ggsave(
  filename = file.path(plot_saving_dir, "all_h_sd.pdf"),
  plot = p_sd_all,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_h_sd.pdf"))

# plot median sd for each sample, color by author
all_h_sd_med <- all_h_sd %>%
  group_by(author, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

((
  p_sd_med <- ggplot(all_h_sd_med, aes(x, y, group = author)) +
    geom_line(aes(color = author), alpha = 0.8, linewidth = 0.3) +
    ylab("Fraction") +
    xlab("Fragment Length (bp)") +
    theme_classic()
  
))

ggsave(
  filename = file.path(plot_saving_dir, "all_h_sd_med.pdf"),
  plot = p_sd_med,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_h_sd_med.pdf"))




# ichor -------------------------------------------------------------------


ichor <- mae[, mae$author == 'V.A. et al', c("sd", "length")]

ichor_long <- longFormat(ichor, colDataCols = c("author", "cohort")) %>%
  tibble::as_tibble() 

ichor_len<- ichor_long %>%
  dplyr::filter(assay == "length") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)

((
p_long <- ggplot(ichor_len, aes(x, y, group = primary)) +
  geom_line(color = "grey", alpha = 0.3, linewidth = 0.2) +
    labs(y = "fraction") +
  theme_classic() +
  facet_wrap(~cohort)
  
  
))


ichor_sd<- ichor_long %>%
  dplyr::filter(assay == "sd") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)

((
  p_sd <- ggplot(ichor_sd, aes(x, y, group = primary)) +
    geom_line(aes(color = cohort), alpha = 0.8, linewidth = 0.2) +
    labs(y = "sd") +
    theme_classic() +
    facet_wrap(~cohort)
  
  
))


p_long / p_sd


# delfi -------------------------------------------------------------------



delfi <- mae[, mae$author == 'S.C. et al', c("sd", "length")]

delfi_long <- longFormat(delfi, colDataCols = "cohort") %>%
  tibble::as_tibble() 

delfi_len<- delfi_long %>%
  dplyr::filter(assay == "length") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)

((
  p_long <- ggplot(delfi_len, aes(x, y, group = primary)) +
    geom_line(color = "grey", alpha = 0.8, linewidth = 0.2) +
    labs(y = "fraction") +
    facet_wrap(~cohort)
  
  
))


delfi_sd<- delfi_long %>%
  dplyr::filter(assay == "sd") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value)

((
  p_sd <- ggplot(delfi_sd, aes(x, y, group = primary)) +
    geom_line(color = "grey", alpha = 0.8, linewidth = 0.2) +
    labs(y = "sd") +
    facet_wrap(~cohort)
  
  
))


p_long / p_sd






# plot sd -----------------------------------------------------------------

sd <- delfi_l_tibble %>%
  dplyr::filter(assay == "sd") %>%
  dplyr::mutate(x = stringr::str_remove(rowname, "^isize_") |> as.numeric()) %>%
  dplyr::mutate(y = value)

((
  
  p_sd <- ggplot(sd, aes(x, y, group = primary )) +
    geom_line(color = "grey", alpha = 0.3, linewidth = 1) +
    facet_wrap(~cohort)
  
))

# plot length -----------------------------------------------------------------

len <- delfi_l_tibble %>%
  dplyr::filter(assay == "length") %>%
  dplyr::mutate(x = stringr::str_remove(rowname, "^isize_") |> as.numeric()) %>%
  dplyr::mutate(y = value)

((
  
  p_len <- ggplot(len, aes(x, y, group = primary )) +
    geom_line(aes(color = cohort), alpha = 0.3, linewidth = 1) +
    facet_wrap(~cohort)
  
))



# plot s/l -----------------------------------------------------------------

sl_ratio <- delfi_l_tibble %>%
  dplyr::filter(assay == "sl_ratio") %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)

((
  
  p_sl_ratio <- ggplot(sl_ratio, aes(x, y, group = primary )) +
    geom_line(color = "grey", alpha = 0.3, linewidth = 1) +
    facet_wrap(~cohort)
  
))


# plot cnv -----------------------------------------------------------------

cnv <- delfi_l_tibble %>%
  dplyr::filter(assay == "cnv") %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)

((
  
  p_sl_cnv <- ggplot(cnv, aes(x, y, group = primary )) +
    geom_line(color = "grey", alpha = 0.3, linewidth = 1) +
    facet_wrap(~cohort)
  
))


# plot motif -----------------------------------------------------------------

motif_ratio <- delfi_l_tibble %>%
  dplyr::filter(assay == "motif_ratio") %>%
  dplyr::mutate(x = rowname) %>%
  dplyr::mutate(y = value)

((
  
  p_motif_ratio <- ggplot(motif_ratio, aes(x, y, group = primary )) +
    geom_line(color = "grey", alpha = 0.3, linewidth = 1) +
    facet_wrap(~cohort) +
    coord_cartesian(ylim = c(0, 3))
  
))


################################################################################
# get all samples
################################################################################

filter_criteria <-  (!mae$author %in% c('M.S. et al', 'Stephen Cristiano', 'K.S. et al')) & (mae$timepoint == 1 | is.na(mae$timepoint))

all <- mae[, filter_criteria, c("sd", "length", "sl_ratio", "cnv", "motif_ratio")]

all_long <- longFormat(all, 
            colDataCols = c("cohort", 
                            "author", 
                            "stage", 
                            "library_kit", 
                            "extraction_kit", 
                            "seq_platform", 
                            "age", 
                            "gender",
                            "timepoint",
                            "clinical_tf",
                            "ichorcna_tf",
                            "mito_frac",
                            "sample_type",
                            "data_source"
                            )) %>%
  tibble::as_tibble() 


all_long_filter <- all_long %>%
  dplyr::filter(timepoint == 1 | is.na(timepoint)) %>%
  dplyr::filter(sample_type == "plasma") 


################################################################################
# length
################################################################################
all_len<- all_long_filter %>%
  dplyr::filter(assay == "length") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value) %>%
  mutate(bicohort = ifelse(cohort == "Healthy", "Healthy", "Cancer"))

# add ichorcna_tf_strat
all_len <- all_len %>%
  mutate(ichorcna_tf_strat = case_when(
    ichorcna_tf <= 0.01 ~ "(0, 0.01]",
    ichorcna_tf <= 0.03 & ichorcna_tf > 0.01 ~ "(0.01, 0.03]",
    ichorcna_tf <= 0.1 & ichorcna_tf > 0.03 ~ "(0.03, 0.1]",
    ichorcna_tf <= 0.2 & ichorcna_tf > 0.1 ~ "(0.1, 0.2]",
    ichorcna_tf <= 1 & ichorcna_tf > 0.2 ~ "(0.2, 1]",
    TRUE ~ "unknown"
  ))


# plot all samples together, color by author
((
  p_long <- ggplot(all_len, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.2, linewidth = 0.2) +
    ylab("Fraction") +
    xlab("Fragment Length (bp)") +
    theme_classic()
  
))

ggsave(
  filename = file.path(plot_saving_dir, "all_len.pdf"),
  plot = p_long,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_len.pdf"))


# plot median length for each bicohor

all_len_med <- all_len %>%
  group_by(bicohort, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

((
  p_long_med <- ggplot(all_len_med, aes(x, y, group = bicohort)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.4) +
    ylab("Median Fraction") +
    xlab("Fragment Length (bp)") +
    theme_classic()
  
))


ggsave(
  filename = file.path(plot_saving_dir, "all_len_med.pdf"),
  plot = p_long_med,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_len_med.pdf"))

# plot length, facet by author
((
  p_long <- ggplot(all_len, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("Fraction") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # set color palette
    scale_color_nord(palette = "aurora") +
    facet_wrap(~author)
))

ggsave(
  filename = file.path(plot_saving_dir, "all_len_facet.pdf"),
  plot = p_long,
  width = 6,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_len_facet.pdf"))



# plot length, facet by author and bicohort
((
  p_long <- ggplot(all_len, aes(x, y, group = primary)) +
    geom_line(aes(color = bicohort), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("Fraction") +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    theme_classic() +
    # hide legend
    theme(legend.position = "none") +
    # set color palette
    scale_color_manual(values = c("Cancer" = "grey", 
                                "Healthy" = "royalblue")) +
    facet_grid( bicohort ~ author)
))

ggsave(
  filename = file.path(plot_saving_dir, "all_len_facet_author_cohort.pdf"),
  plot = p_long,
  width = 10,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_len_facet_author_cohort.pdf"))



# plot length facet by ichorcna_tf_strat
all_len_healthy <- all_len %>%
  filter(bicohort == "Healthy") %>%
  filter(ichorcna_tf_strat != 'unknown') %>%
  mutate(ichorcna_tf_strat = "Healthy")




all_len_cancer <- all_len %>%
  filter(bicohort == "Cancer") %>%
  filter(ichorcna_tf_strat != 'unknown')

# combine the two dataframes
all_len_clean <- rbind(all_len_healthy, all_len_cancer)


label_levels <- c("Healthy", 
                  "(0, 0.01]", 
                  "(0.01, 0.03]", 
                  "(0.03, 0.1]", 
                  "(0.1, 0.2]", 
                  "(0.2, 1]")

# set the order of ichorcna_tf_strat, "Healthy" first
all_len_clean$ichorcna_tf_strat <- factor(all_len_clean$ichorcna_tf_strat, 
                                        levels = label_levels)

# count the number of sample for each ichorcna_tf_strat
make_label <- all_len_clean %>%
  group_by(ichorcna_tf_strat) %>%
  summarize(n = length(unique(primary))) %>%
  ungroup() %>%
  mutate(facet_label = paste0("ichorCNA TF ", ichorcna_tf_strat, "\n", "(n=", n, ")"))

facet_label <- make_label$facet_label
names(facet_label) <- make_label$ichorcna_tf_strat


((
  p_long <- ggplot(all_len_clean, aes(x, y, group = primary)) +
    #geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = 1, fill = "lightgrey", alpha = 0.2)+
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.8, linewidth = 0.2) +
    xlab("Fragment Length (bp)") +
    ylab("Fraction") +
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
  filename = file.path(plot_saving_dir, "all_len_facet_ichorcna_tf_strat.pdf"),
  plot = p_long,
  width = 7,
  height = 4,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "all_len_facet_ichorcna_tf_strat.pdf"))


# median length for each ichorcna_tf_strat

# calculate the median of all samples 
all_len_healthy_med <- all_len_healthy %>%
  group_by(x) %>%
  summarize(y = median(y)) %>%
  ungroup() %>%
  mutate(ichorcna_tf_strat = "Healthy")


all_len_cancer_med <- all_len_cancer %>%
  group_by(ichorcna_tf_strat, x) %>%
  summarize(y = median(y)) %>%
  ungroup()

all_len_med <- rbind(all_len_healthy_med, all_len_cancer_med)

# set the order of ichorcna_tf_strat, "Healthy" first
all_len_med$ichorcna_tf_strat <- factor(all_len_med$ichorcna_tf_strat, 
                                        levels = c("Healthy", 
                                                    "(0, 0.01]", 
                                                    "(0.01, 0.03]", 
                                                    "(0.03, 0.1]", 
                                                    "(0.1, 0.2]", 
                                                    "(0.2, 1]"))

# plot median length for each ichorcna_tf_strat
((
  p_long_med <- ggplot(all_len_med, aes(x, y, group = ichorcna_tf_strat)) +
    geom_rect(xmin = 100, xmax = 155, ymin = 0, ymax = 1, fill = "lightgrey", alpha = 0.2)+
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.95, linewidth = 0.7) +
    geom_vline(xintercept = 167, color = "black", linetype = "dashed", linewidth = 0.1 ) +
    ylab("Median Fraction") +
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
  filename = file.path(plot_saving_dir, "all_len_med_ichorcna_tf_strat.pdf"),
  plot = p_long_med,
  width = 6,
  height = 3,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_len_med_ichorcna_tf_strat.pdf"))


# for each healthy sample, calcualte the total fraction between 100 and 155 bp
all_len_healthy_100_155 <- all_len_healthy %>%
  filter(x >= 100 & x <= 155) %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarize(total_fraction = sum(y)) %>%
  ungroup()

# for each cancer sample, calcualte the total fraction between 100 and 155 bp
all_len_cancer_100_155 <- all_len_cancer %>%
  filter(x >= 100 & x <= 155) %>%
  group_by(primary, ichorcna_tf_strat) %>%
  summarize(total_fraction = sum(y)) %>%
  ungroup()

# combine the two dataframes
all_len_100_155 <- rbind(all_len_healthy_100_155, all_len_cancer_100_155)

# set the order of ichorcna_tf_strat, "Healthy" first
all_len_100_155$ichorcna_tf_strat <- factor(all_len_100_155$ichorcna_tf_strat, 
                                        levels = c("Healthy", 
                                                    "(0, 0.01]", 
                                                    "(0.01, 0.03]", 
                                                    "(0.03, 0.1]", 
                                                    "(0.1, 0.2]", 
                                                    "(0.2, 1]"))

# plot as box plot
((
  p_box <- ggplot(all_len_100_155, aes(x = ichorcna_tf_strat, y = total_fraction)) +
  #boxplot without outliers
    geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA ) +
    geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.55) +
    ylab("Frac of 100~155bp fragments") +
    xlab("ichorCNA Tumor Fraction") +
    # remove legend 
    guides(color = guide_legend(title = NULL)) +
    # increase the legend text size and line size
    theme(legend.text = element_text(size = 8), 
          legend.title = element_text(size = 8),
          legend.key.height = unit(0.9, "cm")) +
    theme_classic() +
    scale_fill_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
))

ggsave(
  filename = file.path(plot_saving_dir, "all_len_100_155_boxplot.pdf"),
  plot = p_box,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(plot_saving_dir, "all_len_100_155_boxplot.pdf"))

# use patchwork to combine the  plots, hide all legends

p_long_med_nolegend <- p_long_med + theme(legend.position = "none")
p_box_nolegend <- p_box + theme(legend.position = "none")

ctdna_is_shorter <- p_long  + 
                    (p_long_med_nolegend / p_box_nolegend) +
                    plot_layout(widths = c(1.3, 1))

ggsave(
  filename = file.path(plot_saving_dir, "ctdna_is_shorter.pdf"),
  plot = ctdna_is_shorter,
  width = 10,
  height = 6,
  dpi = 300
)
message("Saved to ", file.path(plot_saving_dir, "ctdna_is_shorter.pdf"))


################################################################################
# sd 
################################################################################

all_sd <- all_long_filter %>%
  dplyr::filter(assay == "sd") %>%
  mutate(x = stringr::str_remove(rowname, "isize_") |> as.numeric()) %>%
  mutate(y = value) %>%
  mutate(bicohort = ifelse(cohort == "Healthy", "Healthy", "Cancer"))

# add ichorcna_tf_strat
all_sd <- all_sd %>%
  mutate(ichorcna_tf_strat = case_when(
    ichorcna_tf <= 0.01 ~ "(0, 0.01]",
    ichorcna_tf <= 0.03 & ichorcna_tf > 0.01 ~ "(0.01, 0.03]",
    ichorcna_tf <= 0.1 & ichorcna_tf > 0.03 ~ "(0.03, 0.1]",
    ichorcna_tf <= 0.2 & ichorcna_tf > 0.1 ~ "(0.1, 0.2]",
    ichorcna_tf <= 1 & ichorcna_tf > 0.2 ~ "(0.2, 1]",
    TRUE ~ "unknown"
  ))


