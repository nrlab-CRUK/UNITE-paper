library(tidyverse)
library(cfDNAPro)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(rslurm)
library(factoextra)

# print the package version

source("/home/nrlab/wang04/ulyses/plot_mae/gather_outliers.R")
source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")

path <- downsampled_plasma_bams_data_path_hyper
motif_plot_saving_dir_base <- "/home/nrlab/wang04/ulyses/plot_mae/plots/motif_genomewide"
motif_plot_saving_dir <- file.path(motif_plot_saving_dir_base, paste0("iteration_", iteration_num_hyper))
# make dir if not exist
if (!dir.exists(motif_plot_saving_dir)) {
  dir.create(motif_plot_saving_dir, recursive = TRUE)
}

slurm_options_param <- list(partition = "epyc", time = "0-1:00:00", mem = "8G")



read_motif <- function(x) {
  ans <- cfDNAPro::readBam(x, genome_label = "hg19") %>%
    cfDNAPro::callMotif(genome_label = "hg19", motif_type = "s", motif_length = 3L)
  return(ans)
}

read_motif_from_RDS <- function(x) {
  ans <- readRDS(x) %>%
    cfDNAPro::callMotif(genome_label = "hg19", motif_type = "s", motif_length = 3L)
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

bam_id_in_mae <- colData(mae) %>%
  rownames()

fl1 <- list.files(path = path, pattern = pattern1, recursive = TRUE, full.names = TRUE)
fl1_id <- sapply(fl1, get_bam_id)
fl1_filtered <- fl1[which(fl1_id %in% bam_id_in_mae)]
fl2 <- list.files(path = path, pattern = pattern2, recursive = TRUE, full.names = TRUE)
fl2_id <- sapply(fl2, get_bam_id)
fl2_filtered <- fl2[which(fl2_id %in% bam_id_in_mae)]
fl <- c(fl1_filtered, fl2_filtered)
fl_id <- sapply(fl, get_bam_id)
# issue a warning if the number of element in fl_id is different from bam_id_in_mae
if (length(fl_id) != length(bam_id_in_mae)) {
  message("No. of ids in fl_id:", length(fl_id))
  message("No. of ids in mae:", length(bam_id_in_mae))
  # report the difference
  setdiff(fl_id, bam_id_in_mae)
  setdiff(bam_id_in_mae, fl_id)
}


sjob_motif1 <- rslurm::slurm_map(
  x = fl1_filtered |> as.list() |> setNames(fl1_filtered),
  f = read_motif,
  jobname = "rslurm_motif1",
  global_objects = c("read_motif"),
  job_array_task_limit = 150,
  pkgs = c("cfDNAPro", "tidyverse"),
  nodes = 1000,
  cpus_per_node = 1,
  processes_per_node = 1,
  rscript_path = "/home/nrlab/tools/anaconda3/envs/R4_3/bin/Rscript",
  slurm_options = slurm_options_param,
  submit = TRUE
)

fl1_read_motif <- rslurm::get_slurm_out(sjob_motif1, outtype = "raw", wait = TRUE)
ans1 <- bind_rows(fl1_read_motif, .id = "file")


sjob_motif2 <- rslurm::slurm_map(
  x = fl2_filtered |> as.list() |> setNames(fl2_filtered),
  f = read_motif_from_RDS,
  jobname = "rslurm_motif2",
  global_objects = c("read_motif_from_RDS"),
  job_array_task_limit = 150,
  pkgs = c("cfDNAPro", "tidyverse"),
  nodes = 1000,
  cpus_per_node = 1,
  processes_per_node = 1,
  rscript_path = "/home/nrlab/tools/anaconda3/envs/R4_3/bin/Rscript",
  slurm_options = slurm_options_param,
  submit = TRUE
)

fl2_read_motif <- rslurm::get_slurm_out(sjob_motif2, outtype = "raw", wait = TRUE)
ans2 <- bind_rows(fl2_read_motif, .id = "file")


all_motif <- bind_rows(ans1, ans2)


# fl1_read_motif <- bettermc::mclapply(fl1, read_motif, mc.cores = 20)
# fl2_read_motif <- bettermc::mclapply(fl2, read_motif_from_RDS, mc.cores = 20)
# fl_read_motif <- c(fl1_read_motif, fl2_read_motif)
# names(fl_read_motif) <- as.vector(fl)
# all_motif <- bind_rows(fl_read_motif, .id = "file")



all_motif <- all_motif %>%
  mutate(file_dir = dirname(file)) %>%
  mutate(file_name = basename(file))

motif_bam_id <- map(all_motif$file_name, get_bam_id) %>% unlist()
all_motif$bam_id <- motif_bam_id


# save ans object as RDS file under /home/nrlab/wang04/ulyses/plot_mae/plots/motif
saveRDS(all_motif, file.path(motif_plot_saving_dir, "motif_genomewide_s3.rds"))



################################################################################
# sd
################################################################################

# calcualte frac of each motif
mae_col_data <- mae_col_data %>%
  dplyr::mutate(gender = case_when(
    gender == "M" ~ "male",
    gender == "F" ~ "female",
    TRUE ~ gender
  )) %>%
  select(-n)

all_motif_add_meta <- left_join(mae_col_data, all_motif, by = c("primary" = "bam_id"))





################################################################################
# plots
################################################################################


###############################################################################
# pca ONLY healthy
###############################################################################
data_selected_len_cancer <- all_motif_add_meta %>% dplyr::filter(bicohort %in% c("Cancer"))
data_selected_len_healthy <- all_motif_add_meta %>% dplyr::filter(bicohort %in% c("Healthy"))

all_motif_add_meta1 <- select(data_selected_len_healthy, primary, motif, n, author, data_source, bicohort, library_kit, extraction_kit, seq_platform, age, gender)
# expand the "rowname" columns to wider format
all_motif_add_meta2 <- pivot_wider(all_motif_add_meta1, names_from = motif, values_from = n) %>%
  # change primary to rownames
  column_to_rownames("primary")
all_motif_add_meta2_active <- all_motif_add_meta2[, 9:ncol(all_motif_add_meta2)]



pca_plot_title <- paste0("PCA of s3 motifs ")
plot_suffix <- "_pca_motif_healthy.pdf"

groups_bicohort <- all_motif_add_meta2$bicohort
groups_bicohort <- factor(groups_bicohort, levels = c("Healthy", "Cancer"))
groups_author <- all_motif_add_meta2$author
groups_data_source <- all_motif_add_meta2$data_source

groups_library_kit <- all_motif_add_meta2$library_kit
groups_extraction_kit <- all_motif_add_meta2$extraction_kit
groups_seq_platform <- all_motif_add_meta2$seq_platform
groups_age <- all_motif_add_meta2$age
groups_gender <- all_motif_add_meta2$gender

# change NA to character "NA"
groups_library_kit[is.na(groups_library_kit)] <- "Not Available"
groups_extraction_kit[is.na(groups_extraction_kit)] <- "Not Available"
groups_seq_platform[is.na(groups_seq_platform)] <- "Not Available"
groups_age[is.na(groups_age)] <- "Not Available"
groups_gender[is.na(groups_gender)] <- "Not Available"


# make harmonied all_motif_add_meta2_active
all_motif_add_meta2_active_harmonized <- sva::ComBat_seq(all_motif_add_meta2_active |> t(), batch = groups_author) |> t()

# convert to the matrix to fraction of each motif
all_motif_add_meta2_active_frac <- all_motif_add_meta2_active / rowSums(all_motif_add_meta2_active)
all_motif_add_meta2_active_harmonized_frac <- all_motif_add_meta2_active_harmonized / rowSums(all_motif_add_meta2_active_harmonized)

res.pca <- prcomp(all_motif_add_meta2_active_frac)
p_ranking <- fviz_eig(res.pca)

#------------------------------------------------------------------------------

authors_vec <- all_motif_add_meta2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_len_author <- fviz_pca_ind(res.pca,
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
p_pca_len_author <- p_pca_len_author + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_author", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_author,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

#------------------------------------------------------------------------------
# harmonized

pca_plot_title <- paste0("PCA of s3 motifs (harmonized)")
plot_suffix <- "_pca_motif_healthy.pdf"

res.pca.harmonized <- prcomp(all_motif_add_meta2_active_harmonized_frac)
p_ranking_harmonized <- fviz_eig(res.pca.harmonized)


authors_vec <- all_motif_add_meta2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_len_author_harmonized <- fviz_pca_ind(res.pca.harmonized,
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
p_pca_len_author_harmonized <- p_pca_len_author_harmonized + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_author_harmonized", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_author_harmonized,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)



#------------------------------------------------------------------------------
p_pca_len_bicohort <- fviz_pca_ind(res.pca,
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

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_bicohort", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_bicohort,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)



#------------------------------------------------------------------------------

p_pca_len_datasource <- fviz_pca_ind(res.pca,
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

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_datasource", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_datasource,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)


#--------------------------------------------------------------------------
# plot by library_kit
p_pca_len_library_kit <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_library_kit, # color by groups
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Library Kit",
  title = pca_plot_title,
  repel = TRUE
)

shape_values <- if (any(is.na(groups_library_kit))) {
  seq(0, length(unique(groups_library_kit)) + 1)
} else {
  seq(0, length(unique(groups_library_kit)))
}
p_pca_len_library_kit <- p_pca_len_library_kit + scale_shape_manual(values = shape_values)

ggsave(
  filename = file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_library_kit", plot_suffix)),
  plot = p_pca_len_library_kit,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)

message("Saved to ", file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_library_kit", plot_suffix)))

#--------------------------------------------------------------------------
# plot by extraction_kit
p_pca_len_extraction_kit <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_extraction_kit, # color by groups
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Extraction Kit",
  title = pca_plot_title,
  repel = TRUE
)

shape_values <- if (any(is.na(groups_extraction_kit))) {
  seq(0, length(unique(groups_extraction_kit)) + 1)
} else {
  seq(0, length(unique(groups_extraction_kit)))
}

p_pca_len_extraction_kit <- p_pca_len_extraction_kit + scale_shape_manual(values = shape_values)

ggsave(
  filename = file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_extraction_kit", plot_suffix)),
  plot = p_pca_len_extraction_kit,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)

message("Saved to ", file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_extraction_kit", plot_suffix)))


#--------------------------------------------------------------------------
# plot by seq_platform
p_pca_len_seq_platform <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_seq_platform, # color by groups
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Sequencing Platform",
  title = pca_plot_title,
  repel = TRUE
)

p_pca_len_seq_platform <- p_pca_len_seq_platform + scale_shape_manual(values = seq(0, length(unique(groups_seq_platform))))

ggsave(
  filename = file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_seq_platform", plot_suffix)),
  plot = p_pca_len_seq_platform,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)

message("Saved to ", file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_seq_platform", plot_suffix)))

#--------------------------------------------------------------------------

# plot by age
p_pca_len_age <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_age, # color by groups
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Age",
  title = pca_plot_title,
  repel = TRUE
)

ggsave(
  filename = file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_age", plot_suffix)),
  plot = p_pca_len_age,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)

message("Saved to ", file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_age", plot_suffix)))

#--------------------------------------------------------------------------
# plot by gender
p_pca_len_gender <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_gender, # color by groups
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "Age",
  title = pca_plot_title,
  repel = TRUE
)

ggsave(
  filename = file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_gender", plot_suffix)),
  plot = p_pca_len_gender,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)

message("Saved to ", file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_gender", plot_suffix)))






###############################################################################
# pca: healthy + cancer
###############################################################################

all_motif_add_meta1 <- select(all_motif_add_meta, primary, motif, n, author, data_source, bicohort, ichorcna_tf_strat)
#  %>% filter(ichorcna_tf_strat %in% c("Healthy", "(0.2, 1]"))
# expand the "rowname" columns to wider format
all_motif_add_meta2 <- pivot_wider(all_motif_add_meta1, names_from = motif, values_from = n) %>%
  # change primary to rownames
  column_to_rownames("primary")
all_motif_add_meta2_active <- all_motif_add_meta2[, 5:ncol(all_motif_add_meta2)]
# convert to the matrix to fraction of each motif
all_motif_add_meta2_active_frac <- all_motif_add_meta2_active / rowSums(all_motif_add_meta2_active)

pca_plot_title <- paste0("PCA of s3 motifs ")

# plot file suffix
plot_suffix <- "_pca_motif_healthy_and_cancer.pdf"

# fviz_eig(res.pca)
groups_bicohort <- all_motif_add_meta2$bicohort
groups_bicohort <- factor(groups_bicohort, levels = c("Healthy", "Cancer"))
groups_author <- all_motif_add_meta2$author
groups_data_source <- all_motif_add_meta2$data_source
groups_ichorcna_tf_strat <- all_motif_add_meta2$ichorcna_tf_strat



#--------------------------------------------------------------------------
res.pca <- prcomp(all_motif_add_meta2_active_frac)
# plot by ichorcna_tf_strat
p_pca_len_ichorcna_tf_strat <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_ichorcna_tf_strat, # color by groups
  palette = ichorcna_tf_strat_color_hyper[groups_ichorcna_tf_strat |>
    unique() |>
    as.vector()],
  alpha.ind = 0.6,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "ichorCNA TF Strat",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_ichorcna_tf_strat", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_ichorcna_tf_strat,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

# harmonized
pca_plot_title <- paste0("PCA of s3 motifs (Harmonized)")
plot_suffix <- "_pca_motif_healthy_and_cancer.pdf"

all_motif_add_meta2_active_harmonized <- sva::ComBat_seq(all_motif_add_meta2_active |> t(),
  batch = groups_author |> as.vector(),
  group = groups_ichorcna_tf_strat |> as.vector()
) |>
  t()
# convert to fraction
all_motif_add_meta2_active_harmonized_frac <- all_motif_add_meta2_active_harmonized / rowSums(all_motif_add_meta2_active_harmonized)

res.pca_harmonized <- prcomp(all_motif_add_meta2_active_harmonized_frac)
# plot by ichorcna_tf_strat
p_pca_len_ichorcna_tf_strat_harmonized <- fviz_pca_ind(res.pca_harmonized,
  # hide label
  geom = c("point"),
  col.ind = groups_ichorcna_tf_strat, # color by groups
  # palette = ichorcna_tf_strat_color_hyper[groups_ichorcna_tf_strat |> unique() |> as.vector()],
  alpha.ind = 0.8,
  addEllipses = TRUE, # Concentration ellipses
  legend.title = "ichorCNA TF Strat",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_ichorcna_tf_strat_harmonized", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_ichorcna_tf_strat_harmonized,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)




#--------------------------------------------------------------------------

p_pca_len_bicohort <- fviz_pca_ind(res.pca,
  # hide label
  geom = c("point"),
  col.ind = groups_bicohort, # color by groups
  # palette = c("royalblue", "#f95959"),
  palette = c(healthy_color_hyper, cancer_color_hyper),
  alpha.ind = 0.8,
  addEllipses = TRUE, # Concentration ellipses
  ellipse.type = "confidence",
  legend.title = "Cohort",
  title = pca_plot_title,
  repel = TRUE
)

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_bicohort", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_bicohort,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)




#------------------------------------------------------------------------------

authors_vec <- all_motif_add_meta2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_len_author <- fviz_pca_ind(res.pca,
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
p_pca_len_author <- p_pca_len_author + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_author", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_author,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

#------------------------------------------------------------------------------

p_pca_len_datasource <- fviz_pca_ind(res.pca,
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

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_datasource", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_datasource,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)



###############################################################################
# pca: cancer
###############################################################################
data_selected_len_cancer <- all_motif_add_meta %>% dplyr::filter(bicohort %in% c("Cancer"))
data_selected_len_healthy <- all_motif_add_meta %>% dplyr::filter(bicohort %in% c("Healthy"))

all_motif_add_meta1 <- select(data_selected_len_cancer, primary, rowname, value, author, data_source, bicohort, ichorcna_tf_strat)
# expand the "rowname" columns to wider format
all_motif_add_meta2 <- pivot_wider(all_motif_add_meta1, names_from = rowname, values_from = value) %>%
  # change primary to rownames
  column_to_rownames("primary")
all_motif_add_meta2_active <- all_motif_add_meta2[, 5:ncol(all_motif_add_meta2)]
res.pca <- prcomp(all_motif_add_meta2_active)


length_range_label_isize_min <- min(all_motif_add_meta1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label_isize_max <- max(all_motif_add_meta1$rowname |> str_remove("isize_") |> as.numeric())
length_range_label <- paste0("(", length_range_label_isize_min, " - ", length_range_label_isize_max, "bp)")
pca_plot_title <- paste0("PCA of Fragment Length ", length_range_label)

# plot file suffix
plot_suffix <- "_pca_length_cancer.pdf"

fviz_eig(res.pca)
groups_bicohort <- all_motif_add_meta2$bicohort
groups_bicohort <- factor(groups_bicohort, levels = c("Healthy", "Cancer"))
groups_author <- all_motif_add_meta2$author
groups_data_source <- all_motif_add_meta2$data_source
groups_ichorcna_tf_strat <- all_motif_add_meta2$ichorcna_tf_strat

#--------------------------------------------------------------------------
# plot by ichorcna_tf_strat
p_pca_len_ichorcna_tf_strat <- fviz_pca_ind(res.pca,
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

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_ichorcna_tf_strat", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_ichorcna_tf_strat,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)


#--------------------------------------------------------------------------
p_pca_len_bicohort <- fviz_pca_ind(res.pca,
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

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_bicohort", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_bicohort,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)




#------------------------------------------------------------------------------

authors_vec <- all_motif_add_meta2$author |> unique()
authors_color_code <- author_color_values_hyper[authors_vec]

p_pca_len_author <- fviz_pca_ind(res.pca,
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
p_pca_len_author <- p_pca_len_author + scale_shape_manual(values = seq(0, length(author_color_values_hyper)))

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_author", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_author,
  width = 18,
  height = 13,
  dpi = 300,
  units = "cm"
)
message("Saved to ", plotfile)

#------------------------------------------------------------------------------

p_pca_len_datasource <- fviz_pca_ind(res.pca,
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

plotfile <- file.path(motif_plot_saving_dir, paste0("all_motif_add_meta_datasource", plot_suffix))

ggsave(
  filename = plotfile,
  plot = p_pca_len_datasource,
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
# plot CCC motif fract for each ichorcna_tf_strat as boxplot
################################################################################

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
ylab_txt <- "Genomewide CCC/AAA Ratio"
  p_box_motif <- ggplot(data1, aes(x = ichorcna_tf_strat, y = ratio)) +
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
    geom_hline(yintercept = median_healthy_CCC_AAA, color = "black", linetype = "dashed", linewidth = 0.1) +
    ylab(ylab_txt) +
    xlab("ichorCNA Tumor Fraction") +
    # add p-value to the plot
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
    theme_classic() +
  theme(axis.title.x = element_blank()) +
    # remove legend
    theme(legend.position = "none") +
    scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
    # rotate x-axis tick and text by 45 degree
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # axis title and text size to 5
    theme(axis.title = element_text(size = 5.5), axis.text = element_text(size = 5))
  # facet_wrap(~author)



ggsave(
  filename = file.path(motif_plot_saving_dir, paste("boxplot_motif_CCCAAA_ratio.pdf", sep = "")),
  plot = p_box_motif,
  width = 4.27,
  height = 4.72,
  units = "cm",
  dpi = 300
)

message("Saved to ", file.path(motif_plot_saving_dir, paste("boxplot_motif_CCCAAA_ratio.pdf", sep = "")))

################################################################################
# CGN_NCG
################################################################################

# plot as box plot
ylab_txt <- "Genomewide CGN/NCG Ratio"
p_box_motif <- ggplot(CGN_NCG, aes(x = ichorcna_tf_strat, y = CGN_NCG_ratio)) +
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
  geom_hline(yintercept = median_healthy_CGN_NCG, color = "black", linetype = "dashed", linewidth = 0.1) +
  ylab(ylab_txt) +
  xlab("ichorCNA Tumor Fraction") +
  theme_classic() +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", ref.group = "Healthy", vjust = 0.5) +
  theme_classic() +
  # remove x axis title
  theme(axis.title.x = element_blank()) +
  # remove legend
  theme(legend.position = "none") +
  scale_fill_manual(values = ichorcna_tf_strat_color_hyper[ichorcna_tf_strat_hyper]) +
  # remove legend
  # scale_fill_nord(palette = "aurora") +
  # rotate x-axis tick and text by 45 degree
  theme(axis.title = element_text(size = 5.5), axis.text = element_text(size = 5))

p_box_motif_facet_author <- p_box_motif + facet_wrap(~author)
p_box_motif_facet_cohort <- p_box_motif + facet_wrap(~cohort)




ggsave(
  filename = file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio.pdf", sep = "")),
  plot = p_box_motif,
  width = 4.27,
  height = 4.72,
  units = "cm",
  dpi = 300
)
message("Saved to ", file.path(motif_plot_saving_dir, paste("boxplot_CGN_NCG_ratio.pdf", sep = "")))

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

cleanup_files(sjob_motif1, wait = TRUE)
cleanup_files(sjob_motif2, wait = TRUE)
