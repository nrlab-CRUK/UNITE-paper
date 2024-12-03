
# Install pacman if not already installed
if (!require("pacman")) install.packages("pacman")

# Load multiple packages
pacman::p_load(tidyverse, 
patchwork, 
MultiAssayExperiment,
nord,
plotly,
htmlwidgets,
ggpubr)

message("Loading MAE RDS data")
source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
mae <- readRDS(mae_file_latest)
replace_xgb_input_data <- FALSE
replace_test_set <- FALSE
################################################################################
# get all samples
################################################################################

# prepare filtering criteria
author_filter <- mae$author %in% authors_hyper
timepoint_filter <- !is.na(mae$timepoint) & mae$timepoint == 1
sample_type_filter <- mae$sample_type == sample_type_hyper
ichorcna_tf_filter <- !is.na(mae$ichorcna_tf)
cohort_filter <- mae$cohort %in% cohorts_hyper
outlier_filter <- !rownames(colData(mae)) %in% outliers_hyper

filter_criteria <- author_filter &
  timepoint_filter &
  sample_type_filter &
  ichorcna_tf_filter &
  cohort_filter &
  outlier_filter

# figure out which assay to keep
assay_to_keep <- assay_to_keep_hyper

# perform filtering
all <- mae[, filter_criteria, assay_to_keep]

# conver to long formats
all_long <- longFormat(all,
  colDataCols = c(
    "cohort",
    "patient_id",
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
  )
) %>%
  tibble::as_tibble()

all_long_filter <- all_long %>%
  mutate(bicohort = case_when(
    cohort == "Healthy" ~ "Healthy",
    cohort %in% noncancer_disease_hyper ~ "Non-cancer Disease",
    TRUE ~ "Cancer"
  )) %>%
  filter(!is.na(ichorcna_tf)) %>%
  filter(bicohort != "Non-cancer Disease")


# add ichorcna_tf_strat
message("Adding ichorcna_tf_strat")
all_long_filter_h <- all_long_filter %>%
  filter(bicohort == "Healthy") %>%
  mutate(ichorcna_tf_strat = "Healthy")

all_long_filter_c <- all_long_filter %>%
  filter(bicohort == "Cancer")
# mutate(ichorcna_tf_strat = case_when(
#   ichorcna_tf <= 0.03 ~ "[0, 0.03]",
#   # ichorcna_tf <= 0.03 & ichorcna_tf > 0.01 ~ "(0.01, 0.03]",
#   ichorcna_tf <= 0.1 & ichorcna_tf > 0.03 ~ "(0.03, 0.1]",
#   ichorcna_tf <= 0.2 & ichorcna_tf > 0.1 ~ "(0.1, 0.2]",
#   ichorcna_tf <= 1 & ichorcna_tf > 0.2 ~ "(0.2, 1]",
#   TRUE ~ "unknown"
# ))

all_long_filter_c <- add_ichorcna_tf_strat_col(all_long_filter_c)


all_long_filter <- rbind(all_long_filter_h, all_long_filter_c)

# set the factor levels of ichorcna_tf_strat
all_long_filter$ichorcna_tf_strat <- factor(all_long_filter$ichorcna_tf_strat,
  levels = ichorcna_tf_strat_col_levels_hyper
)


################################################################################
# tidy up the meta data for samples
################################################################################

# stages--------------------------------
all_long_filter$stage <- all_long_filter$stage %>%
  dplyr::case_match(
    "0" ~ "Others",
    "X" ~ "Others",
    "IA" ~ "I",
    "IB" ~ "I",
    "IIA" ~ "II",
    "IIB" ~ "II",
    "IIIA" ~ "III",
    "IIIB" ~ "III",
    "IIIC" ~ "III",
    "IVA" ~ "IV",
    "IVB" ~ "IV",
    .default = all_long_filter$stage
  )

# set factor order to stage
all_long_filter$stage <- factor(all_long_filter$stage, levels = names(stage_color_hyper))

# library_kit--------------------------------
all_long_filter$library_kit <- all_long_filter$library_kit %>%
  dplyr::case_match(
    "Kapa Hyper Prep kit" ~ "KAPA Hyper Prep Kit",
    "KAPA Hyper Prep kit" ~ "KAPA Hyper Prep Kit",
    "Kapa Library Preparation Kit" ~ "KAPA Library Preparation Kit",
    "Thruplex DNA−seq" ~ "ThruPLEX DNA−seq",
    "Thruplex DNA−seq " ~ "ThruPLEX DNA−seq",
    "Thruplex Plasma−seq" ~ "ThruPLEX Plasma−seq",
    "Thruplex Plasma−seq " ~ "ThruPLEX Plasma−seq",
    .default = all_long_filter$library_kit
  )


# extraction_kit--------------------------------


all_long_filter$extraction_kit <- all_long_filter$extraction_kit %>%
  dplyr::case_match(
    "Qiagen Circulating DNA kit" ~ "Qiagen Circulating DNA Kit",
    "QIAsymphony DSP Circulating DNA kit" ~ "QIAsymphony DSP Circulating DNA Kit",
    .default = all_long_filter$extraction_kit
  )
# seq_platform--------------------------------


all_long_filter$seq_platform <- all_long_filter$seq_platform %>%
  dplyr::case_match(
    "Illumina NextSeq or MiSeq" ~ "NextSeq or MiSeq",
    .default = all_long_filter$seq_platform
  )

# Author--------------------------------

all_long_filter$author <- all_long_filter$author %>%
  dplyr::case_match(
    "A.S. et al" ~ "Santonja et al, 2023",
    "A.Z. et al" ~ "Zviran et al, 2020",
    "P.U. et al" ~ "Ulz et al, 2019",
    "F.M. et al" ~ "Mouliere et al, 2018",
    "K.H. et al" ~ "Gale et al, 2022",
    "K.S. et al" ~ "Sun et al, 2019",
    "M.S. et al" ~ "Snyder et al, 2016",
    "P.J. et al" ~ "Jiang et al, 2015",
    "P.M. et al" ~ "This study",
    "P.P. et al" ~ "Peneder et al, 2021",
    "S.C. et al" ~ "Cristiano et al, 2019",
    "V.A. et al" ~ "Adalsteinsson et al, 2017",
    .default = all_long_filter$author
  )

author_levels <- c(
  "Santonja et al, 2023",
  "Zviran et al, 2020",
  "Ulz et al, 2019",
  "Mouliere et al, 2018",
  "Gale et al, 2022",
  "Sun et al, 2019",
  "Snyder et al, 2016",
  "Jiang et al, 2015",
  "This study",
  "Peneder et al, 2021",
  "Cristiano et al, 2019",
  "Adalsteinsson et al, 2017"
)
# set factor order to author
all_long_filter$author <- factor(all_long_filter$author, levels = author_levels)


# save
if (replace_xgb_input_data) {
  saveRDS(object = all_long_filter, file = xgboost_input_data_file_hyper)
  message("Saved to: ", xgboost_input_data_file_hyper)
}


# tidy up mae_col_data

mae_col_data <- all_long_filter %>%
  # remove 'patient_id'
  dplyr::select(-patient_id) %>%
  group_by(primary, colname, cohort, author, stage, library_kit, extraction_kit, seq_platform, age, gender, timepoint, clinical_tf, ichorcna_tf, mito_frac, sample_type, data_source, bicohort, ichorcna_tf_strat) %>%
  summarise(n = n()) %>%
  ungroup()





################################################################################
# make final test set
################################################################################

if (replace_test_set) {
  set.seed(seed_number_hyper)
  final_test <- all_long_filter %>%
    select(-rowname, -value, -assay) %>%
    # remove duplicated rows
    distinct() %>%
    group_by(author, cohort, ichorcna_tf_strat) %>%
    slice_sample(replace = FALSE, prop = final_test_prop_hyper)

  # count how many samples are in the final test set, grouped by author, cohort,  ichorcna_tf_strat and assay
  final_test_count <- final_test %>%
    group_by(author, cohort, ichorcna_tf_strat) %>%
    summarise(n = n()) %>%
    ungroup()

  final_test_primary_vec <- final_test[["primary"]]
  final_test_patient_id_vec <- final_test[["patient_id"]]
  save(final_test,
    final_test_count,
    final_test_primary_vec,
    final_test_patient_id_vec,
    file = final_test_hyper
  )
  message("Saved to: ", final_test_hyper)
}
