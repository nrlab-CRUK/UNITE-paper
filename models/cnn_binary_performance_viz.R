library(tidyverse)
library(nord)
library(argparse)
library(gridExtra)
library(grid)
library(kableExtra)

# write wd as a parameter for this file

parser <- argparse::ArgumentParser(description = "Params parser")
# add parameter 'wd', set default as current working directory
parser$add_argument('--wd', type = 'character', help = 'working directory', default = getwd())
# add paramter 'file_pattern', set default as "performance.*\\.debug.csv$"
parser$add_argument('--file_pattern', type = 'character', help = 'file pattern', default = "performance.*\\.debug.csv$")
# add paramter 'metrics_plot_file', set default as "ulyses_binary_performance_TF_stratification.pdf"  
parser$add_argument('--metrics_plot_file', type = 'character', help = 'metrics plot file', default = "ulyses_binary_performance_TF_stratification.pdf")
# add parameter 'metrics_to_plot', set default as c("accuracy", "AUROC", "sen_95spe", "sen_98spe", "sen_99spe")
parser$add_argument('--metrics_to_plot', type = 'character', help = 'metrics to plot', default = c("accuracy", "AUROC", "sen_95spe", "sen_98spe", "sen_99spe"))
# add ichorcna_tf_strat as a parameter, set default as c( "(0, 0.01]", "(0.01, 0.03]", "(0.03, 0.1]", "(0.1, 0.2]", "(0.2, 1]", "all")
parser$add_argument('--ichorcna_tf_strat', type = 'character', help = 'ichorcna_tf_strat', default = c( "(0, 0.01]", "(0.01, 0.03]", "(0.03, 0.1]", "(0.1, 0.2]", "(0.2, 1]", "all"))

# add plot_width as a parameter, set default as 10
parser$add_argument('--plot_width', type = 'numeric', help = 'plot width', default = 10)
# add plot_height as a parameter, set default as 5
parser$add_argument('--plot_height', type = 'numeric', help = 'plot height', default = 5)
# add plot_units as a parameter, set default as "in"
parser$add_argument('--plot_units', type = 'character', help = 'plot units', default = "in")
# add plot_dpi as a parameter, set default as 300
parser$add_argument('--plot_dpi', type = 'numeric', help = 'plot dpi', default = 300)


args <- parser$parse_args()
wd <- args$wd
file_pattern <- args$file_pattern
ichorcna_tf_strat <- args$ichorcna_tf_strat

plot_width <- args$plot_width
plot_height <- args$plot_height
plot_units <- args$plot_units
plot_dpi <- args$plot_dpi

metrics_plot_file_fullname <- file.path(wd, args$metrics_plot_file)
metrics_to_plot <- args$metrics_to_plot
table_pdf_file <- file.path(wd, "performance_table.pdf")


fl <- list.files(path = wd, pattern = file_pattern, full.names = TRUE)

fl_read <- map(fl, read_csv, show_col_types = FALSE) %>% 
  map(pivot_longer, cols = everything(), names_to = "metrics", values_to = "value") %>%
  setNames(fl) %>%
  bind_rows(.id = "file") %>%
  mutate(`ichorcna_tf_strat` = case_when(
    str_detect(file, "0_0.01") ~ "(0, 0.01]",
    str_detect(file, "0.01_0.03") ~ "(0.01, 0.03]",
    str_detect(file, "0.03_0.1") ~ "(0.03, 0.1]",
    str_detect(file, "0.1_0.2") ~ "(0.1, 0.2]",
    str_detect(file, "0.2_1") ~ "(0.2, 1]",
    str_detect(file, "0_1") ~ "all",
  ))

# set order of ichorcna_strat
fl_read$`ichorcna_tf_strat` <- factor(fl_read$`ichorcna_tf_strat`, levels = ichorcna_tf_strat)




fl_read_train_test <- fl_read %>% 
  mutate(train_test = case_when(
    str_detect(metrics, "train") ~ "train",
    str_detect(metrics, "test") ~ "test",
    TRUE ~ "NA"
  )) %>%
  mutate(metrics = str_remove(metrics, "^train_|^test_"))

fl_read_train_test$train_test <- factor(fl_read_train_test$train_test, levels = c("train", "test"))

# filter tibble using metric_to_plot
fl_read_train_test_use <- fl_read_train_test %>% filter(metrics %in% metrics_to_plot)

# save fl_read_train_test_use as rds file
fl_read_train_test_use_rds_file <- file.path(wd, "performance_tibble.rds")
saveRDS(fl_read_train_test_use, fl_read_train_test_use_rds_file)
message("Tibble Saved to ", fl_read_train_test_use_rds_file)

#box plot comparing two datasets

p <- ggplot(fl_read_train_test_use, aes(x = `ichorcna_tf_strat`, y = value)) +
  # add boxplot, do not show outliers
  geom_boxplot(aes(fill = `ichorcna_tf_strat`), outlier.shape = NA, color = "#2a2727") +
  # add median value as text annotation beside the boxplot, align the text to y = 1
  #stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 3)), hjust = 0) +
  geom_jitter(color = "black", fill = "#e9e6e6", alpha = 0.3, width = 0.2, size = 0.3) +
  labs(x = 'ichorCNA Tumor Fraction', y = "Performance") +
  theme_classic() +
  # hide legend
  theme(legend.position = "none") +
  # increase the legend text size and line size
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8),
        legend.key.height = unit(0.9, "cm")) +
  scale_fill_nord(palette = "aurora") +
  # rotate x-axis tick and text by 45 degree
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # facet by metrics and train_test
  facet_grid(train_test ~ metrics, scales = "free_y") +
  # set facet label background color to light grey
  theme(strip.background = element_rect(fill = "#e9e6e6", color = "#e9e6e6")) 


# save plot
ggsave(
  metrics_plot_file_fullname,
  plot = p,
  width = plot_width,
  height = plot_height,
  units = plot_units,
  dpi = plot_dpi)

message("Saved to ", metrics_plot_file_fullname)



# summarize as table

summary_tibble <- fl_read_train_test_use %>% 
  group_by(`ichorcna_tf_strat`, train_test, metrics) %>%
  summarize(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(`to_report` = paste(round(mean, 3), round(sd, 3), sep = "±")) %>%
  select(-c(mean, sd)) %>% 
  pivot_wider(names_from = `ichorcna_tf_strat` , values_from = to_report) %>%
  # rename 'train_test' as ''
  rename(` ` = train_test) 


table_obj <- kable(summary_tibble, "html") %>%
  kable_styling()

# save table as png
save_kable(table_obj, table_pdf_file)
message("Saved to ", table_pdf_file)




