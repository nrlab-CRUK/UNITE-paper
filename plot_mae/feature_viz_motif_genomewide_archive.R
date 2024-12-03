library(tidyverse)
library(cfDNAPro)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)

path <- "/mnt/scratchc/nrlab/wang04/ulyses/data"
motif_plot_saving_dir <- "/home/nrlab/wang04/ulyses/plot_mae/plots/motif"

source("/home/nrlab/wang04/ulyses/plot_mae/gather_outliers.R")


read_motif <- function(x) {
  ans <- cfDNAPro::readBam(x, genome_label = "hg19") %>%
    cfDNAPro::callMotif(genome_label = "hg19", motif_type = "s", motif_length = 4L)
  return(ans)
}

read_motif_from_RDS <- function(x) {
  ans <- readRDS(x) %>%
    cfDNAPro::callMotif(genome_label = "hg19", motif_type = "s", motif_length = 4L)
  return(ans)
}

get_bam_id <- function(x) {
  if (stringr::str_detect(x, "EE\\d+\\.hg")) {
    bam_id <- stringr::str_extract(basename(x), "(^EE\\d+)\\.hg\\d\\d", group = 1)
  } else {
    bam_id <- stringr::str_extract(basename(x), "(^(SLX|DL)\\-\\d+\\.(\\w+\\d+(\\w+)?(\\-\\w+\\d+(\\w+)?)?)|(\\w+))\\.", group = 1)
  }
  return(bam_id)
}
# pattern

# bamfile pattern
pattern1 <- "\\.0\\.1x\\.mrkdup\\.bam$"
# finaleDB pattern
pattern2 <- "hg19\\.frag\\.tsv\\.bgz\\.GRanges\\.rds\\.1M\\.rds$"

fl1 <- list.files(path = path, pattern = pattern1, recursive = TRUE, full.names = TRUE)
fl2 <- list.files(path = path, pattern = pattern2, recursive = TRUE, full.names = TRUE)
fl <- c(fl1, fl2)

fl1_read_motif <- bettermc::mclapply(fl1, read_motif, mc.cores = 20)
fl2_read_motif <- bettermc::mclapply(fl2, read_motif_from_RDS, mc.cores = 20)

fl_read_motif <- c(fl1_read_motif, fl2_read_motif)

names(fl_read_motif) <- as.vector(fl)

all_motif <- bind_rows(fl_read_motif, .id = "file")



all_motif <- all_motif %>%
  mutate(file_dir = dirname(file)) %>%
  mutate(file_name = basename(file))

motif_bam_id <- map(all_motif$file_name, get_bam_id) %>% unlist()
all_motif$bam_id <- motif_bam_id


# save ans object as RDS file under /home/nrlab/wang04/ulyses/plot_mae/plots/motif
saveRDS(all_motif, file.path(plot_saving_dir, "motif_s4.rds"))


################################################################################
# get all samples
################################################################################
source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")

################################################################################
# sd
################################################################################

# calcualte frac of each motif

all_motif <- all_motif %>%
  group_by(bam_id) %>%
  mutate(total_frag = sum(n)) %>%
  ungroup() %>%
  mutate(frac = n / total_frag) %>%
  filter(!str_detect(file, "csf"))

all_motif_add_meta <- left_join(mae_col_data, all_motif, by = c("primary" = "bam_id"))




################################################################################
# plots
################################################################################



# plot CCC motif fract for each ichorcna_tf_strat as boxplot

data1 <- all_motif_add_meta %>%
  select(primary, ichorcna_tf_strat, motif, frac, author, cohort) %>%
  filter(motif %in% c("CCC", "AAA")) %>%
  pivot_wider(names_from = motif, values_from = frac) %>%
  mutate(ratio = CCC / AAA)

# calculate the median of all samples
median_healthy_CCC_AAA <- median(data1$ratio[data1$ichorcna_tf_strat == "Healthy"])

# NCG motif

NCG <- all_motif_add_meta %>%
  select(primary, ichorcna_tf_strat, motif, frac, author, cohort) %>%
  filter(str_detect(motif, "[ACGT]CG")) %>%
  group_by(primary, ichorcna_tf_strat, author, cohort) %>%
  summarize(NCG = sum(frac))

CGN <- all_motif_add_meta %>%
  select(primary, ichorcna_tf_strat, motif, frac, author, cohort) %>%
  filter(str_detect(motif, "CG[ATCG]")) %>%
  group_by(primary, ichorcna_tf_strat, author, cohort) %>%
  summarize(CGN = sum(frac))

CGN_NCG <- inner_join(NCG, CGN, by = c("primary", "ichorcna_tf_strat", "author", "cohort")) %>%
  mutate(CGN_NCG_ratio = CGN / NCG)

median_healthy_CGN_NCG <- median(CGN_NCG$CGN_NCG_ratio[CGN_NCG$ichorcna_tf_strat == "Healthy"])
# plot as box plot
ylab_txt <- "CCC/AAA Ratio"
((
  p_box_motif <- ggplot(data1, aes(x = ichorcna_tf_strat, y = ratio)) +
    geom_violin(fill = "lightgrey", color ="lightgrey", width = 1.1) +
    # boxplot without outliers
    geom_jitter(color = "black", 
      fill = "#e9e6e6", 
      alpha = 0.9, 
      width = 0.2, 
      size = 0.2, 
      shape = 21,
      stroke = 0.01) +
    geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA, width = 0.2) +
    geom_hline(yintercept = median_healthy_CCC_AAA, color = "black", linetype = "dashed", linewidth = 0.1) +
    ylab(ylab_txt) +
    xlab("ichorCNA Tumor Fraction") +
    # add p-value to the plot
    stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
    theme_classic() +
    # remove legend
    theme(legend.position = "none") +
    scale_fill_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # axis title and text size to 5
    theme(axis.title = element_text(size = 6), axis.text = element_text(size = 5.5))
  # facet_wrap(~author)

))


ggsave(
  filename = file.path(motif_plot_saving_dir, paste("boxplot_motif.pdf", sep = "")),
  plot = p_box_motif,
  width = 6.5,
  height = 5,
  units = "cm",
  dpi = 300
)

message("Saved to ", file.path(motif_plot_saving_dir, paste("boxplot_motif.pdf", sep = "")))




# plot as box plot
ylab_txt <- "5'-CGN/NCG Ratio"
p_box_motif <- ggplot(CGN_NCG, aes(x = ichorcna_tf_strat, y = CGN_NCG_ratio)) +
geom_violin(fill = "lightgrey", color ="lightgrey", alpha = 0.5) +
  # boxplot without outliers
  geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.55) +
  geom_hline(yintercept = median_healthy_CGN_NCG, color = "black", linetype = "dashed", linewidth = 0.1) +
  geom_boxplot(aes(fill = ichorcna_tf_strat), outlier.shape = NA) +
  ylab(ylab_txt) +
  xlab("ichorCNA Tumor Fraction") +
  theme_classic() +
  # remove legend
  theme(legend.position = "none") +
  scale_fill_nord(palette = "aurora") +
  # rotate x-axis tick and text by 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # add p-value to the plot
  stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy")

p_box_motif_facet_author <- p_box_motif + facet_wrap(~author)
p_box_motif_facet_cohort <- p_box_motif + facet_wrap(~cohort)




ggsave(
  filename = file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio.pdf", sep = "")),
  plot = p_box_motif,
  width = 6,
  height = 4,
  dpi = 300
)
ggsave(
  filename = file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio_facet_author.pdf", sep = "")),
  plot = p_box_motif_facet_author,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  filename = file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio_facet_cohort.pdf", sep = "")),
  plot = p_box_motif_facet_cohort,
  width = 10,
  height = 10,
  dpi = 300
)

message("Saved to ", file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio.pdf", sep = "")))

message("Saved to ", file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio_facet_author.pdf", sep = "")))

message("Saved to ", file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio_facet_cohort.pdf", sep = "")))




################################################################################
# plot motif frac as line plot, facet by ichorcna_tf_strat
################################################################################


# plot as line plot
ylab_txt <- "Fraction"
((
  p_line_motif <- ggplot(all_motif_add_meta, aes(x = motif, y = frac, group = primary)) +
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.8, linewidth = 0.2) +
    ylab(ylab_txt) +
    xlab("Fragment Length (bp)") +
    # add p-value to the plot
    # stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy") +
    theme_classic() +
    # remove legend
    theme(legend.position = "none") +
    scale_color_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # x-axis text smaller
    theme(axis.text.x = element_text(size = 3)) +
    facet_wrap(~ichorcna_tf_strat)

))

ggsave(
  filename = file.path(motif_plot_saving_dir, paste("lineplot_motif.pdf", sep = "")),
  plot = p_line_motif,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(motif_plot_saving_dir, paste("lineplot_motif.pdf", sep = "")))

################################################################################
# plot motif median frac as line plot, facet by ichorcna_tf_strat, with median
################################################################################
# calculate the median of all samples
all_motif_med <- all_motif_add_meta %>%
  group_by(ichorcna_tf_strat, motif) %>%
  summarize(median_frac = mean(frac)) %>%
  ungroup()

# plot as line plot, color by ichorcna_tf_strat

ylab_txt <- "Median Fraction"
((
  p_line_motif_med <- ggplot(all_motif_med, aes(x = motif, y = median_frac, group = ichorcna_tf_strat)) +
    geom_line(aes(color = ichorcna_tf_strat), alpha = 0.98, linewidth = 0.4) +
    ylab(ylab_txt) +
    xlab("5' motif (3-mer)") +
    theme_classic() +
    # remove legend
    # theme(legend.position = "none") +
    # remove legend title
    theme(legend.title = element_blank()) +
    scale_color_nord(palette = "aurora") +
    # rotate x-axis tick and text by 45 degree
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # x-axis text smaller
    theme(axis.text.x = element_text(size = 3))

))

ggsave(
  filename = file.path(motif_plot_saving_dir, paste("lineplot_motif_median.pdf", sep = "")),
  plot = p_line_motif_med,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved to ", file.path(motif_plot_saving_dir, paste("lineplot_motif_median.pdf", sep = "")))
