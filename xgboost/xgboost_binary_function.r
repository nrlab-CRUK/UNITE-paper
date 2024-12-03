
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidymodels))
library(rslurm)
library(ggpubr)
library(nord)
library(patchwork)
library(vip)
library(argparse)

###############################################################################
# parameters
###############################################################################

# use argparse to parse the parameters
parser <- argparse::ArgumentParser()
# add xgboost_input_data_file argument
parser$add_argument('--xgboost_input_data_file', type = 'character', help = 'xgboost_input_data_file', default = "xgboost_data_all.rda")
# add outdir argument, default as current directory
parser$add_argument('--outdir', type = 'character', help = 'outdir', default = './')
# add n_fold argument, default as 5
parser$add_argument('--n_fold', type = 'integer', help = 'n_fold', default = 5)
# add n_repeat argument, default as 1
parser$add_argument('--n_repeat', type = 'integer', help = 'n_repeat', default = 10)
# add partition argument, default as "general"
parser$add_argument('--partition', type = 'character', help = 'partition', default = "general")

# get the parameters
args <- parser$parse_args()
xgboost_input_data_file <- args$xgboost_input_data_file
outdir <- args$outdir
n_fold <- args$n_fold
n_repeat <- args$n_repeat
partition <- args$partition

# get host name and based on that set the partition 

rslurm_options <- list(time = '15:00:00', mem = '32G' , partition = partition)

# load the data
load(xgboost_input_data_file)

csv_file_folder <- paste(outdir, tf_label, sep = '')

  # create if not exist
  if (!dir.exists(csv_file_folder)) {
  dir.create(csv_file_folder)
}


###############################################################################

my_model <- function(feat_label, feat_list = feat_list, n_fold = 5, n_repeat = 20){

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidymodels))



p_tune_file <- paste("/home/nrlab/wang04/ulyses/xgboost/", feat_label, "/xgb_tune.png", sep = '')
p_vip_file <- paste("/home/nrlab/wang04/ulyses/xgboost/", feat_label, "/xgb_vip.png", sep = '')
p_roc_file <- paste("/home/nrlab/wang04/ulyses/xgboost/", feat_label, "/xgb_roc.png", sep = '')


vb_split <- initial_split(feat_list[[feat_label]], prop = 0.70, strata = bicohort)
vb_train <- training(vb_split)
vb_test <- testing(vb_split)


xgb_spec <- boost_tree(
  trees = 1000,
  tree_depth = tune(), min_n = tune(),
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune()                          ## step size
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")


xgb_grid <- grid_latin_hypercube(
  tree_depth(),
  min_n(),
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), vb_train),
  learn_rate(),
  size = 30
)

xgb_wf <- workflow() %>%
  add_formula(bicohort ~ .) %>%
  add_model(xgb_spec)


vb_folds <- vfold_cv(data = vb_train, v = n_fold,  repeats = n_repeat, strata = bicohort)
doParallel::registerDoParallel()

set.seed(234)
message("start tuning")
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)

best_auc <- select_best(xgb_res, "roc_auc")

final_xgb <- finalize_workflow(
  xgb_wf,
  best_auc
)
message("last fit")
final_res <- last_fit(final_xgb, vb_split)
test_metrics <- collect_metrics(final_res)
test_metrics_list <- list(test_metrics)
names(test_metrics_list) <- feat_label



###############################################################################
# plots
###############################################################################
message("plotting")
p_tune<- xgb_res %>%
  collect_metrics() %>%
  filter(.metric == "roc_auc") %>%
  select(mean, mtry:sample_size) %>%
  pivot_longer(mtry:sample_size,
               values_to = "value",
               names_to = "parameter"
  ) %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC") +
  theme_classic()

ggsave(p_tune_file, p_tune, width = 10, height = 8, units = "in")



p_vip <- final_xgb %>%
  fit(data = vb_train) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")
p_vip <- p_vip + theme_classic()

ggsave(p_vip_file, p_vip, width = 10, height = 8, units = "in")



#plot roc
p_roc <- final_res %>%
  collect_predictions() %>%
  roc_curve(bicohort, .pred_Cancer, event_level = "second") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(linewidth = 1.5, color = "midnightblue") +
  geom_abline(
    lty = 2, alpha = 0.5,
    color = "gray50",
    linewidth = 1.2
  )+ 
  theme_classic()
ggsave(p_roc_file, p_roc, width = 10, height = 8, units = "in")

message("return test metrics")
return(test_metrics_list)
message("done")




	
}

# for loop in R

fl <- list()
for (i in 1:50){

  repeat_name <- paste("repeat", i, sep = '')
  message(paste(repeat_name, "start..."))
  plot_file_name <- file.path(csv_file_folder, paste(repeat_name, "_tf", tf_label, ".pdf", sep = ''))
  csv_file_name <- file.path(csv_file_folder, paste(repeat_name, "_tf", tf_label, ".csv", sep = ''))
  job <- slurm_map(x=  label_list, f = my_model, feat_list = feat_list, n_fold = 5, n_repeat = 20, 
          jobname = paste(repeat_name, "_tf", tf_label, sep = ''),
          nodes = 1000,
          cpus_per_node = 1,
          processes_per_node = 1,
          slurm_options = rslurm_options,
          #global_objects = global_objects,
          submit = TRUE
          )




results <- get_slurm_out(job, outtype = 'raw', wait = TRUE)


ans <- flatten(results) %>%
  bind_rows(.id = "feat") %>%
  group_by(.metric) %>%
  arrange(.estimate, .by_group = TRUE) %>%
  # factor feat
  mutate(feat = factor(feat, levels = feat))
fl[i] <- ans
write_csv(ans, csv_file_name)

# ggplot facet by .metric

p_final <- ggplot(ans, aes(x = feat, y = .estimate, color = .metric )) +
  geom_point() +
  facet_wrap(~.metric ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Feature", y = "Performance")

ggsave(plot_file_name, 
  p_final,
  width = 8, height = 4)

# clean up job files
cleanup_files(job)


}
