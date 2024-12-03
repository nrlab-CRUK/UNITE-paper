
library(ggupset)
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(argparse)
library(kableExtra)
library(gridExtra)
library(grid)
library(nord)



# add parameters
parser <- argparse::ArgumentParser()
# add wd argument, default as current directory
parser$add_argument('--wd', type = 'character', help = 'wd', default = ".")

# parse args
args <- parser$parse_args()
wd <- args$wd

plot_file <- file.path(wd, "final_boxplot.pdf")
table_pdf_file <- file.path(wd, "summary_table.pdf")
###############################################################################
# plot the performance
###############################################################################

fl <- list.files(path = wd,
                 "repeat\\d+_.*\\.csv", 
		 recursive = TRUE,
                 full.names=TRUE)
fl_read <- lapply(fl, read_csv, show_col_types = FALSE)
names(fl_read) <- fl
fl_read_bind <- bind_rows(fl_read, .id = "id") %>%
  select(-c(.estimator, .config)) %>%
  mutate(ichorcna_strat = case_when(
	stringr::str_detect(id, "0-0.01") ~ '(0, 0.01]',
	stringr::str_detect(id, "0.01-0.03") ~ '(0.01, 0.03]',
	stringr::str_detect(id, "0.03-0.1") ~ '(0.03, 0.1]',
	stringr::str_detect(id, "0.1-0.2") ~ '(0.1, 0.2]',
	stringr::str_detect(id, "0.2-1") ~ '(0.2, 1]',
	stringr::str_detect(id, "all") ~ 'all',
	TRUE ~ NA_character_
  ))
# set order of ichorcna_strat
fl_read_bind$ichorcna_strat <- factor(fl_read_bind$ichorcna_strat, levels = c('all', '(0, 0.01]', '(0.01, 0.03]', '(0.03, 0.1]', '(0.1, 0.2]', '(0.2, 1]'))


# plot the feature performance of 'cnv_sd_length_motif'
ans <- fl_read_bind %>%
filter(feat == 'cnv_sd_length_motif') %>%
select(-c(id, feat))


# summarize as table

summary_tibble <- ans %>% 
  group_by(ichorcna_strat, .metric) %>%
  summarize(mean = mean(.estimate), sd = sd(.estimate), .groups = "drop") %>%
  mutate(`to_report` = paste(round(mean, 3), round(sd, 3), sep = "±")) %>%
  select(-c(mean, sd)) %>% 
  pivot_wider(names_from = `ichorcna_strat`, values_from = to_report) %>%
  rename(metrics = .metric) %>%
  # change 'roc_auc' to 'AUROC'
  mutate(metrics = ifelse(metrics == 'roc_auc', 'AUROC', metrics))


table_obj <- kable(summary_tibble, "html") %>%
  kable_styling()

# save table as png
save_kable(table_obj, table_pdf_file)
message(table_pdf_file)


