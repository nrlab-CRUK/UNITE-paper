


library(tidyverse)

wd <- "/archive/Groups/NRLab/fs07/Shared_BamFiles/sWGS/LUCID_plasma_KH/collapsed_bams"
output_dir <- "/home/nrlab/wang04/ulyses/meta_data/lucid"
lucid_meta_file <- file.path(output_dir, "lucid_meta.csv")

healthy_dir <- file.path(wd, "healthy_control")
cancer_dir <- wd

healthy_bc_file <- file.path(wd, "healthy_BC.txt")

healthy_bc <- read_tsv(healthy_bc_file, col_names = FALSE) %>%
  pull(1)
healthy_file_list <- list.files(healthy_dir, pattern = "bam$", full.names = TRUE)
cancer_file_list <- list.files(cancer_dir, pattern = "bam$", full.names = TRUE)
cancer_bc <- str_extract(basename(cancer_file_list), "^SLX-\\d+_SXTHS.\\d+|SLX-\\d+_D\\d{3}_D\\d{3}") %>%
  unique() 


healthy_tibble <- tibble(bam_id = healthy_bc, cohort = "Healthy")
cancer_tidble <- tibble(bam_id = cancer_bc, cohort = "Lung")

lucid_tibble <- bind_rows(healthy_tibble, cancer_tidble) %>%
  mutate(patient_id = paste("LUCID", row_number(), sep = "_")) %>%
  # sustitute the _ to . in the bam_id column
  mutate(bam_id = str_replace(bam_id, "_", ".")) %>%
  mutate(bam_id = str_replace(bam_id, "_", "-"))


write_csv(lucid_tibble, lucid_meta_file)
