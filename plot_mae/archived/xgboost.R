library(tidyverse)
library(tidymodels)
library(rslurm)

source("/home/nrlab/wang04/ulyses/plot_mae/feature_viz_pre_filtering.R")




###############################################################################
# sd 
###############################################################################

all_sd <- all_long_filter %>%
  dplyr::filter(assay == "sd") %>%
  dplyr::select(bicohort, ichorcna_tf,ichorcna_tf_strat, primary, rowname, value)

sd_input <- all_sd %>%
  dplyr::filter(rowname %in% paste("isize_", seq(100, 170), sep = '')) %>%
  pivot_wider(names_from = rowname, values_from = value)

###############################################################################
# length
###############################################################################

all_length <- all_long_filter %>%
  dplyr::filter(assay == "length") %>%
  dplyr::select(bicohort, ichorcna_tf,ichorcna_tf_strat, primary, rowname, value)

length_input <- all_length %>%
  dplyr::filter(rowname %in% paste("isize_", seq(100, 170), sep = '')) %>%
  pivot_wider(names_from = rowname, values_from = value)
###############################################################################
# sl_ratio
###############################################################################

all_sl_ratio <- all_long_filter %>%
  dplyr::filter(assay == "sl_ratio") %>%
  dplyr::select(bicohort, ichorcna_tf,ichorcna_tf_strat, primary, rowname, value)

sl_ratio_input <- all_sl_ratio %>%
  dplyr::filter(rowname != "6p_6" ) %>%
  pivot_wider(names_from = rowname, values_from = value)

###############################################################################
# cnv
###############################################################################

all_cnv <- all_long_filter %>%
  dplyr::filter(assay == "cnv") %>%
  dplyr::select(bicohort, ichorcna_tf,ichorcna_tf_strat, primary, rowname, value)

cnv_input <- all_cnv %>%
  dplyr::filter(rowname != "6p_6" ) %>%
  pivot_wider(names_from = rowname, values_from = value)



###############################################################################
# motif_ratio
###############################################################################

all_motif_ratio <- all_long_filter %>%
  dplyr::filter(assay == "motif_ratio") %>%
  dplyr::select(bicohort, ichorcna_tf,ichorcna_tf_strat, primary, rowname, value)

motif_ratio_input <- all_motif_ratio %>%
  dplyr::filter(rowname != "6p_6" ) %>%
  pivot_wider(names_from = rowname, values_from = value)


###############################################################################
# build xgboost model
###############################################################################

prepare_clean_input <- function(input, tf_strat_vec) {
  input_clean <- input %>%
    dplyr::filter(ichorcna_tf_strat %in% tf_strat_vec) %>%
    dplyr::select( -c(primary,ichorcna_tf_strat, ichorcna_tf))
  
  input_clean$bicohort <- factor(input_clean$bicohort, levels = c("Healthy", "Cancer"))
  
  return(input_clean)
}

tf_strat_vec <- c("Healthy",
                  "(0, 0.01]",
                  "(0.01, 0.03]",
                  "(0.03, 0.1]",
                  "(0.1, 0.2]")

cnv_input_clean <- prepare_clean_input(cnv_input, tf_strat_vec)
sd_input_clean <- prepare_clean_input(sd_input, tf_strat_vec)
length_input_clean <- prepare_clean_input(length_input, tf_strat_vec)
sl_ratio_input_clean <- prepare_clean_input(sl_ratio_input, tf_strat_vec)
motif_ratio_input_clean <- prepare_clean_input(motif_ratio_input, tf_strat_vec)


# model

input_clean <- sl_ratio_input_clean
input_clean <- sd_input_clean
input_clean <- cnv_input_clean

input_clean <- length_input_clean

input_clean <- motif_ratio_input_clean

# set seed
set.seed(123)
vb_split <- initial_split(input_clean, strata = bicohort)
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


vb_folds <- vfold_cv(data = vb_train, v = 5,  repeats = 1, strata = bicohort)

# tune
doParallel::registerDoParallel()

set.seed(234)
xgb_res <- tune_grid(
  xgb_wf,
  resamples = vb_folds,
  grid = xgb_grid,
  control = control_grid(save_pred = TRUE)
)



collect_metrics(xgb_res)


p <- xgb_res %>%
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
  labs(x = NULL, y = "AUC")

ggsave("xgboost/length/xgb_tune.png", p, width = 10, height = 8, units = "in")

best_auc <- select_best(xgb_res, "roc_auc")

final_xgb <- finalize_workflow(
  xgb_wf,
  best_auc
)


library(vip)

p_vip <- final_xgb %>%
  fit(data = vb_train) %>%
  extract_fit_parsnip() %>%
  vip(geom = "point")

ggsave("xgboost/length/xgb_vip.png", p_vip, width = 10, height = 8, units = "in")

final_res <- last_fit(final_xgb, vb_split)

collect_metrics(final_res)


#plot roc
p_roc <- final_res %>%
  collect_predictions() %>%
  roc_curve(bicohort, .pred_Healthy) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(size = 1.5, color = "midnightblue") +
  geom_abline(
    lty = 2, alpha = 0.5,
    color = "gray50",
    size = 1.2
  )
ggsave("xgboost/length/xgb_roc.png", p_roc, width = 10, height = 8, units = "in")
