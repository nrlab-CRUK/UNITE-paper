library(DataEditR)
library(tidyverse)
library(kableExtra)
library(openxlsx)
library(knitr)
library(magrittr)

ulyses_base_path <-"/home/nrlab/wang04/ulyses"
setwd("/home/nrlab/wang04/ulyses")
source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

latex_table_env <- "Stable"
table_save_path_manuscript <- "/home/nrlab/wang04/maintext-Ulysses_nat_biotechnol_maintext"
table_save_path_thesis <- "/home/nrlab/wang04/cruk_ci_phd_thesis/Tables"

supplementary_tables_excel_file <- file.path(table_save_path_manuscript, "Supplementary tables.xlsx")







################################################################################
# handling stable_data_selection
################################################################################

latex_manuscript_outfile <- file.path(table_save_path_manuscript, "plasma_sample_summary_table_train_test_all.tex")
latex_thesis_outfile <- file.path(table_save_path_thesis, "plasma_sample_summary_table_train_test_all.tex")

csv_manuscript_outfile <- file.path(table_save_path_manuscript, "plasma_sample_summary_table_train_test_all.csv")
csv_thesis_outfile <- file.path(table_save_path_thesis, "plasma_sample_summary_table_train_test_all.csv")

# kept data id
latest_metafile <- "/home/nrlab/wang04/ulyses/meta_data/colData_batch1_2.csv"
selected_data <- readRDS("/home/nrlab/wang04/ulyses/xgboost/xgboost_model_data.rds")
selected_bamid <- unique(selected_data$primary)
# test data id
load("/home/nrlab/wang04/ulyses/models/final_test_samples.rda", test_data_rda_env <- new.env())
test_bamid <- test_data_rda_env$final_test_primary_vec




# # copy the ulyses colData csv to phd thesis working directory
# file.copy(from = latest_metafile, to = work_dir, overwrite = TRUE)

# readin the colData csv file
ulyses_colData <- read_csv(latest_metafile, show_col_types = FALSE)

# summarize - plasma the table for thesis
colData_summary_plasma <- ulyses_colData %>%
  filter(sample_type == "plasma") %>%
  filter(author != "Stephen Cristiano") %>%
  filter(author %in% names(author_color_values_hyper)) %>%
  mutate(selected = ifelse(bam_id %in% selected_bamid, 1, 0)) %>%
  mutate(test = ifelse(bam_id %in% test_bamid, 1, 0)) %>%
  group_by(cohort, author) %>%
  summarise(
    n = n(),
    n_selected = sum(selected),
    n_test = sum(test)
  ) %>%
  ungroup() %>%
  rename("Cohort" = "cohort") %>%
  mutate(n_train = n_selected - n_test)

# rename the author names


colData_summary_plasma$author <- author_rename_vec_hyper[colData_summary_plasma$author]
# set the factor levels
colData_summary_plasma$author <- factor(colData_summary_plasma$author, levels = author_rename_vec_hyper)

# stop if NA exist in author column
if (any(is.na(colData_summary_plasma$author))) {
  stop("NA exist in author column")
}

colData_summary_plasma2 <- colData_summary_plasma %>%
  # mutate(frac_selected = n_selected/n) %>%
  mutate(percent_selected = scales::percent(n_selected / n, accuracy = 0.1)) %>%
  mutate(cell_label = paste0(n_train, "+", n_test, "/", n, " (", percent_selected, ")")) %>%
  group_by(Cohort) %>%
  mutate(cohort_n = sum(n)) %>%
  mutate(cohort_n_test = sum(n_test)) %>%
  mutate(cohort_n_train = sum(n_train)) %>%
  mutate(cohort_n_selected = sum(n_selected)) %>%
  mutate(cohort_frac_selected = cohort_n_selected / cohort_n) %>%
  mutate(cohort_percent_selected = scales::percent(cohort_frac_selected, accuracy = 0.1)) %>%
  mutate(cohort_label = paste0(Cohort, ", ", cohort_n_train, "+", cohort_n_test, "/", cohort_n, " (", cohort_percent_selected, ")")) %>%
  group_by(author) %>%
  mutate(author_n = sum(n)) %>%
  mutate(author_n_train = sum(n_train)) %>%
  mutate(author_n_test = sum(n_test)) %>%
  mutate(author_n_selected = sum(n_selected)) %>%
  mutate(author_frac_selected = author_n_selected / author_n) %>%
  mutate(author_percent_selected = scales::percent(author_frac_selected, accuracy = 0.1)) %>%
  mutate(author_label = paste0(author, ", ", author_n_train, "+", author_n_test, "/", author_n, " (", author_percent_selected, ")")) %>%
  ungroup() %>%
  mutate(Cohort = factor(Cohort)) %>%
  mutate(Cohort = fct_relevel(
    Cohort,
    c("Healthy", "Cirrhosis", "Hepatitis B", "Inflammatory Bowel Disease", "Liver Transplant", "Systemic Lupus Erythematosus")
  )) %>%
  arrange(Cohort)

# set the order of the cohort_label
colData_summary_plasma2$cohort_label <- factor(colData_summary_plasma2$cohort_label,
  levels = unique(colData_summary_plasma2$cohort_label)
)


author_label_order <- colData_summary_plasma2 %>%
  arrange(author) %>%
  pull(author_label) %>%
  unique()

colData_summary_plasma3 <- colData_summary_plasma2 %>%
  select(cell_label, cohort_label, author_label)


# pivot wider, column names are author names -----------------------------------
stable_data_selection <- colData_summary_plasma3 %>%
  arrange(cohort_label) %>%
  pivot_wider(
    names_from = author_label,
    values_from = cell_label,
    values_fill = list(cell_label = "")
  ) %>%
  rename("Cohort" = cohort_label) %>%
  select(c("Cohort", all_of(author_label_order))) %T>%
  write_csv(csv_manuscript_outfile) %>%
  write_csv(csv_thesis_outfile)

# start to make the latex table ------------------------------------------------
n_total <- sum(colData_summary_plasma$n)
n_total_train <- sum(colData_summary_plasma$n_train)
n_total_test <- sum(colData_summary_plasma$n_test)
n_total_selected <- sum(colData_summary_plasma$n_selected)
n_total_selected_percentage <- scales::percent(n_total_selected / n_total, accuracy = 0.1)
caption_label <- paste0(n_total_train, "+", n_total_test, "/", n_total, " (", n_total_selected_percentage, ")")
caption <- paste0(
  "Plasma samples selected for model building\n",
  "selected train+selected test/all (Percentage)"
)

# export latex
tmp <- kbl(stable_data_selection,
  booktabs = T, format = "latex",
  label = "plasma_dataset",
  caption = caption,
  table.envir = latex_table_env 
) %>%
  pack_rows("Non-cancer Disease", 2, 6) %>%
  pack_rows("Cancer", 7, nrow(stable_data_selection)) %>%
  add_header_above(c(" " = 1, "NRLAB" = 4, "EGA" = 3, "FINALEDB" = 5)) %>%
  kable_styling(
     latex_options = c("striped", "hold_position", "scale_down")
     ) %>%
  landscape() 
#%T>%
#readr::write_lines(file = latex_manuscript_outfile) %>%
# readr::write_lines(file = latex_thesis_outfile)




################################################################################
# handling stable4
################################################################################








source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
s1_file_data <-  file.path(table_save_path_manuscript, "supplementary_table_s1.csv")
s2_file_data <-  meta_csv_latest
s3_file_data <-  file.path(table_save_path_manuscript, "plasma_sample_summary_table_train_test_all.csv")
s4_file_data <-  file.path(table_save_path_manuscript, "supplementary_table_s4.csv")
s5_file_data <-  file.path(ulyses_base_path, "explore_packaged", "resources", "bins_anno_tiled.rds")
s6_file_data <-  file.path(table_save_path_manuscript, "supplementary_table_s6.csv")
s7_file_data <-  file.path("/scratchc/nrlab/wang04/ulyses/data/finaledb/hg19_downsampled_rds/EE87922.hg19.frag.tsv.bgz.GRanges.rds.1M.rds.ulyses_bind_olap_chunks.rds")
s8_file_data <-  file.path(ulyses_base_path, "0_data_split", "data_split_cohort_split_table.csv")
s9_file_data <-  file.path(table_save_path_manuscript, "supplementary_table_s9.csv")
s10_file_data <-  file.path(table_save_path_manuscript, "supplementary_table_s10.csv")

s1_read <- read_csv(s1_file_data)
s2_read <- read_csv(s2_file_data)
s3_read <- read_csv(s3_file_data)
s4_read <- read_csv(s4_file_data)
s5_read <- readRDS(s5_file_data) |>  pluck(1)
s6_read <- read_csv(s6_file_data)
s7_read <- readRDS(s7_file_data) |>pluck(1) |>slice(1:500)
s8_read <- read_csv(s8_file_data)
s9_read <- read_csv(s9_file_data)
s10_read <- read_csv(s10_file_data)


s1_cap_label <-  "supplementary_table_meta_variable"
s2_cap_label <- "supplementary_table_all_meta"
s3_cap_label <- "supplementary_table_samples"
s4_cap_label <- "supplementary_table_bin_anno"
s5_cap_label <- "supplementary_table_bin_anno_example"
s6_cap_label <- "supplementary_table_frag_anno"
s7_cap_label <- "supplementary_table_bin-frag_overlap"
s8_cap_label <- "supplementary data splitting"
s9_cap_label <- "supplementary software"
s10_cap_label <- "supplementary deep learning software"


s1_cap_full <- "\\textbf{Meta-variables for reach sample.}"
s1_cap_short <- "Meta-variables for reach sample"

s2_cap_full <- "\\textbf{Meta information of all samples analysed in this study.} This table is provided as a separate worksheet in the “\\textit{Supplementary tables.xlsx}” file.” file."
s2_cap_short <- "Meta information of samples."

s3_cap_full <- "\\textbf{The number of plasma samples selected for each cohort from each of the studies.} The number of training, testing and all samples are shown. This table is provided as a separate worksheet in the “\\textit{Supplementary tables.xlsx}” file."
s3_cap_short <- "The number of plasma samples selected for each cohort from each of the studies."

s4_cap_full <- "\\testbf{Bin annotation variables.}"
s4_cap_short <- "Bin annotation variables"


s5_cap_full <- "\\textbf{An example for 100kb bin annotations.} This table is provided as a separate worksheet in the “\\textit{Supplementary tables.xlsx}” file."
s5_cap_short <- "An example for 100kb annotations."

s6_cap_full <- "\\textbf{Fragment annotation variables.}"
s6_cap_short <- "Fragment annotation variables."


s7_cap_full <- "\\textbf{Intersection of bin and fragment annotations.} Due to the large number of rows in the table (i.e., around 1 million fragments per sample with depth of 0.1x ), a truncated form is provided as a separate worksheet in the “\\textit{Supplementary tables.xlsx}” file for simplicity purposes."
s7_cap_short <- "Intersection of bin and fragment annotations."


s8_cap_full <- "\\textbf{The number of samples of each cohort in training and testing data split.} The models were trained independently for each ichorCNA TF stratum."
s8_cap_short <- "Number of samples from each cohort in data split."

s9_cap_full <- "\\textbf{Software, packages, modules and tools used in the data pre-processing.} The upper part of the table shows the software adopted by the Trim Align Pipeline (\\href{https://github.com/nrlab-CRUK/TAP}{TAP}) built for sequencing data trimming and alignment by NRLAB. The lower part clarifies the R packages used in analysis steps such as bin and fragment filtering, annotation and visualisation."

s9_cap_short <- "Software, packages, modules and tools used in data pre-processing."


s10_cap_full <- "\\textbf{Software and Python modules used in model training.}"
s10_cap_short <- "Software and Python modules used in model training."


s1_latex <- kbl(s1_read,
		format = "latex",
		caption = s1_cap_full,
		caption.short = s1_cap_short,
		label = s1_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% kable_styling(latex_options = c("striped", "hold_position"))

#s2_latex <- kbl(s2_read,
#		format = "latex",
#		caption = s2_cap_full,
#		caption.short = s1_cap_short,
#		label = s2_cap_label,
#  		table.envir = latex_table_env,
#		booktabs = TRUE, linesep = "") %>% kable_styling(latex_options = c("striped", "hold_position"))

s3_latex <- kbl(s3_read,
		format = "latex",
		caption = s3_cap_full,
		caption.short = s3_cap_short,
		label = s3_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% kable_styling(latex_options = c("striped", "hold_position"))

s4_latex <- kbl(s4_read,
		format = "latex",
		caption = s4_cap_full,
		caption.short = s4_cap_short,
		label = s4_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% kable_styling(latex_options = c("striped", "hold_position"))

#s5_latex <- kbl(s5_read,
#		format = "latex",
#		caption = s5_cap_full,
#		caption.short = s5_cap_short,
#		label = s5_cap_label,
#  		table.envir = latex_table_env,
#		booktabs = TRUE, linesep = "") %>% kable_styling(latex_options = c("striped", "hold_position"))

s6_latex <- kbl(s6_read,
		format = "latex",
		caption = s6_cap_full,
		caption.short = s6_cap_short,
		label = s6_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% 
    kable_styling(latex_options = c("striped", "hold_position"))

s7_latex <- kbl(s7_read,
  format = "latex",
  caption = s7_cap_full,
  caption.short = s7_cap_short,
  label = s7_cap_label,
  table.envir = latex_table_env,
  booktabs = TRUE, linesep = ""
) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

s8_latex <- kbl(s8_read,
		format = "latex",
		caption = s8_cap_full,
		caption.short = s8_cap_short,
		label = s8_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% 
    kable_styling(latex_options = c("striped", "hold_position", "scale_down")) %>%
    landscape()

s9_latex <- kbl(s9_read,
		format = "latex",
		caption = s9_cap_full,
		caption.short = s9_cap_short,
		label = s9_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% 
    kable_styling(latex_options = c("striped", "hold_position")) %>%
    column_spec(1, width = "6cm") %>%
    column_spec(2, width = "2cm") %>%
    column_spec(3, width = "6cm")

s10_latex <- kbl(s10_read,
		format = "latex",
		caption = s10_cap_full,
		caption.short = s10_cap_short,
		label = s10_cap_label,
  		table.envir = latex_table_env,
		booktabs = TRUE, linesep = "") %>% 
    kable_styling(latex_options = c("striped", "hold_position")) %>%
    column_spec(1, width = "2cm") %>%
    column_spec(2, width = "2cm") %>%
    column_spec(3, width = "8cm")


readr::write_lines(s1_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S1.tex"))
readr::write_lines(s1_latex, file = file.path(table_save_path_thesis,"Supplementary Table S1.tex"))
#readr::write_lines(s2_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S2.tex"))
#readr::write_lines(s2_latex, file = file.path(table_save_path_thesis,"Supplementary Table S2.tex"))
readr::write_lines(s3_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S3.tex"))
readr::write_lines(s3_latex, file = file.path(table_save_path_thesis,"Supplementary Table S3.tex"))
readr::write_lines(s4_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S4.tex"))
readr::write_lines(s4_latex, file = file.path(table_save_path_thesis,"Supplementary Table S4.tex"))
#readr::write_lines(s5_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S5.tex"))
#readr::write_lines(s5_latex, file = file.path(table_save_path_thesis,"Supplementary Table S5.tex"))
readr::write_lines(s6_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S6.tex"))
readr::write_lines(s6_latex, file = file.path(table_save_path_thesis,"Supplementary Table S6.tex"))
readr::write_lines(s7_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S7.tex"))
readr::write_lines(s7_latex, file = file.path(table_save_path_thesis,"Supplementary Table S7.tex"))
readr::write_lines(s8_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S8.tex"))
readr::write_lines(s8_latex, file = file.path(table_save_path_thesis,"Supplementary Table S8.tex"))
readr::write_lines(s9_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S9.tex"))
readr::write_lines(s9_latex, file = file.path(table_save_path_thesis,"Supplementary Table S9.tex"))
readr::write_lines(s10_latex, file = file.path(table_save_path_manuscript,"Supplementary Table S10.tex"))
readr::write_lines(s10_latex, file = file.path(table_save_path_thesis,"Supplementary Table S10.tex"))




table_list <- list(
		   s1 = s1_read,
		   s2 = s2_read,
		   s3 = s3_read,
		   s4 = s4_read,
		   s5 = s5_read,
		   s6 = s6_read,
		   s7 = s7_read,
		   s8 = s8_read,
		   s9 = s9_read,
		   s10 = s10_read)



write.xlsx(x = table_list, file = supplementary_tables_excel_file)
message("Data written to ", supplementary_tables_excel_file)





################################################################################
# commit and push changes
################################################################################

setwd(table_save_path_manuscript)
system("git pull --rebase ; git add *; git commit -m 'update thesis and manuscript tables';git push")

setwd(table_save_path_thesis)
system("git pull --rebase ; git add *; git commit -m 'update thesis and manuscript tables';git push")
################################################################################

