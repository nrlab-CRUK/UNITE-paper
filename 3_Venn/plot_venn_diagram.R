library(ggupset)
library(ggplot2)
library(tidyverse, warn.conflicts = FALSE)
library(argparse)
library(ggVennDiagram)
###############################################################################
# add parameters
###############################################################################

# add parameters
parser <- argparse::ArgumentParser()
# add wd_cnn argument, default as current directory
parser$add_argument("--wd_xgb", type = "character", help = "wd_xgb", default = "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb10")
parser$add_argument("--wd_cnn", type = "character", help = "wd_cnn", default = "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn14")
parser$add_argument("--xgb_model_chosen", type = "character", help = "xgb_model_chosen", default = "xgboost")
parser$add_argument("--xgb_all_feat_name", type = "character", help = "xgb_all_feat_name", default = "cnv_sd_length_ctRatio_slRatio")
parser$add_argument("--cnn_model_chosen", type = "character", help = "cnn_model_chosen", default = "resnet18")
parser$add_argument("--layer_combination_chosen", type = "character", help = "layer_combination_chosen", default = "C4")
# add manual_colors, default as c("#8fbcbb", "#d088c0", "#e0869a", "#eed197", "#8f8cd4")
parser$add_argument("--manual_colors", type = "character", help = "manual_colors", default = c("#8fbcbb", "#d088c0", "#e0869a", "#eed093", "#8f8cd4"))
# add plot_width, default as 200
parser$add_argument("--plot_width", type = "numeric", help = "plot_width", default = 120)
# add plot_height, default as 100
parser$add_argument("--plot_height", type = "numeric", help = "plot_height", default = 50)
# add plot_unit, default as "mm"
parser$add_argument("--plot_unit", type = "character", help = "plot_unit", default = "mm")
# add labs_x, default as "ichorCNA TF Strata"
parser$add_argument("--labs_x", type = "character", help = "labs_x", default = NULL)
# add labs_y, default as "Mean Performance"
parser$add_argument("--labs_y", type = "character", help = "labs_y", default = "Mean performance and 95% CI")

# add all_hight_color, set default as "#b5dae6"
parser$add_argument("--all_hight_color", type = "character", help = "all_hight_color", default = "#c1dde7")



# parse args
args <- parser$parse_args()
wd_cnn <- args$wd_cnn
wd_xgb <- args$wd_xgb
manual_colors <- args$manual_colors
plot_width <- args$plot_width
plot_height <- args$plot_height
plot_unit <- args$plot_unit
labs_x <- args$labs_x
labs_y <- args$labs_y
all_hight_color <- args$all_hight_color

cnn_model_chosen <- args$cnn_model_chosen
layer_combination_chosen <- args$layer_combination_chosen
xgb_model_chosen <- args$xgb_model_chosen
xgb_all_feat_name <- args$xgb_all_feat_name


plot_path <- "/home/nrlab/wang04/ulyses/3_Venn"
xgb_raw_csv <- "/home/nrlab/wang04/ulyses/3_Venn/xgb_raw.csv"
cnn_raw_csv <- "/home/nrlab/wang04/ulyses/3_Venn/cnn_raw.csv"


source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")



###############################################################################
# parse the cnn data
###############################################################################

# get indtest scores

fl_indtest <- list.files(
        path = wd_cnn,
        "repeat_0_fold_0_metrics.csv", recursive = TRUE, full.names = TRUE
)

fl_indtest_read <- bettermc::mclapply(fl_indtest, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_indtest_read) <- fl_indtest

fl_indtest_read_bind <- bind_rows(fl_indtest_read, .id = "id") %>%
        mutate(ichorcna_strat = case_when(
                stringr::str_detect(id, "0_0\\.03") ~ "[0, 0.03]",
                stringr::str_detect(id, "0\\.03_0\\.1") ~ "(0.03, 0.1]",
                stringr::str_detect(id, "0\\.1_1") ~ "(0.1, 1]",
                stringr::str_detect(id, "0_1") ~ "all",
                TRUE ~ "unknown"
        ))
# set order of ichorcna_strat
fl_indtest_read_bind$ichorcna_strat <- factor(fl_indtest_read_bind$ichorcna_strat,
        levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all")
)

# make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
fl_indtest_read_bind <- fl_indtest_read_bind %>%
        mutate(feat = stringr::str_extract(id, "(?<=_model_)C\\d+(_mismatch)?"))

# filter out feat contains "mismatch"
fl_indtest_read_bind <- fl_indtest_read_bind %>%
        filter(!str_detect(feat, "mismatch"))


# plot the feature performance of 'cnv_sd_length_motif'
fl_indtest_read_bind <- fl_indtest_read_bind %>%
        # filter out ichorcna_strat == "unknown"
        filter(ichorcna_strat != "unknown") %>%
        # test_sen_98spe is 1 when test_auroc is 1
        # mutate(test_sen_98spe = ifelse(test_auroc == 1, 1, test_sen_98spe)) %>%
        # pivot_longer to make .metric as a column, all cols start with 'test_' to ".metric", values to ".estimate"
        pivot_longer(cols = starts_with("test_"), names_to = ".metric", values_to = ".estimate")

fl_indtest_read_bind_all_feature_model <- fl_indtest_read_bind %>%
        filter(feat == layer_combination_chosen) %>%
        filter(model_name == cnn_model_chosen)

test_99spec_thresholds <- fl_indtest_read_bind_all_feature_model %>%
        filter(.metric == "test_threshold_99spe") %>%
        select(ichorcna_strat, .metric, .estimate) %>%
        distinct()

test_95spec_thresholds <- fl_indtest_read_bind_all_feature_model %>%
        filter(.metric == "test_threshold_95spe") %>%
        select(ichorcna_strat, .metric, .estimate) %>%
        distinct()

thresholds <- fl_indtest_read_bind_all_feature_model %>%
        select(ichorcna_strat, .metric, .estimate) %>%
        distinct()

###############################################################################

meta_data <- read_csv(meta_csv_latest)

fl <- list.files(
        path = wd_cnn,
        "repeat_0_fold_0_y_test_y_pred_z_test.csv", recursive = TRUE, full.names = TRUE
)

fl_read <- bettermc::mclapply(fl, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_read) <- fl

fl_read_bind <- bind_rows(fl_read, .id = "id") %>%
        mutate(ichorcna_strat = case_when(
                stringr::str_detect(id, "0_0\\.03") ~ "[0, 0.03]",
                stringr::str_detect(id, "0\\.03_0\\.1") ~ "(0.03, 0.1]",
                stringr::str_detect(id, "0\\.1_1") ~ "(0.1, 1]",
                stringr::str_detect(id, "0_1") ~ "all",
                TRUE ~ "unknown"
        ))
# set order of ichorcna_strat
fl_read_bind$ichorcna_strat <- factor(fl_read_bind$ichorcna_strat,
        levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all")
)

# make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
fl_read_bind <- fl_read_bind %>%
        mutate(feat = stringr::str_extract(id, "(?<=_model_)C\\d+(_mismatch)?"))

# filter out feat contains "mismatch"
fl_read_bind <- fl_read_bind %>%
        filter(!str_detect(feat, "mismatch"))


# save fl_read_bind to a csv
fl_read_bind %>%
        write_csv(file.path(plot_path, "cnn_fl_scores_read_bind_indtest.csv"))


# only keep the all_feature model, by str_detect(id, "cnv_sd_length_ctRatio_slRatio")
fl_read_bind <- fl_read_bind %>%
        filter(str_detect(feat, layer_combination_chosen)) %>%
        filter(model_name == cnn_model_chosen)


# add meta_data$cohort to ans2, based on the bam_id
meta_data_selected <- meta_data %>%
        select(bam_id, cohort)

ans2 <- left_join(fl_read_bind, meta_data_selected, by = c("bam_id_test" = "bam_id"))

# y_test as factor
ans2$y_test <- factor(ans2$y_test, levels = c(0, 1))

# recode the y_test to "Healthy" and "Cancer"
ans2$y_test <- recode(ans2$y_test, "0" = "Healthy", "1" = "Cancer")

# make the cohort as factor
ans2$cohort <- factor(ans2$cohort)

# set the level of cohort, make "Healthy" the first level using relevel
ans2$cohort <- relevel(ans2$cohort, ref = "Healthy")

# report the number of samples in each ichorcna_strat and cohort
ans2_summary_cnn <- ans2 %>%
        group_by(ichorcna_strat, cohort) %>%
        summarize(n = n())

ans2 %>%
        group_by(ichorcna_strat, y_test) %>%
        summarize(n = n())

df <- expand.grid(ichorcna_strat_param = unique(ans2$ichorcna_strat))


ichorcna_strat_param_vec <- df$ichorcna_strat_param
# ans2


annotate_detected <- function(ichorcna_strat_param,
                              thresholds_param,
                              ans2_param,
                              cnn_model_chosen_param) {
        # plot
        target_ichorcna_strat <- ichorcna_strat_param
        target_99threshold <- thresholds %>%
                filter(.metric == "test_threshold_99spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)

        target_95threshold <- thresholds %>%
                filter(.metric == "test_threshold_95spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)

        sen_99spec <- thresholds %>%
                filter(.metric == "test_sen_99spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)

        sen_99_spec_label <- sen_99spec %>%
                round(., 4) %>%
                # convert to percentage and append %
                scales::percent(., accuracy = 0.01) %>%
                paste0("Sensitivity = ", .)
        sen_95spec <- thresholds %>%
                filter(.metric == "test_sen_95spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)
        sen_95_spec_label <- sen_95spec %>%
                round(., 4) %>%
                # convert to percentage and append %
                scales::percent(., accuracy = 0.01) %>%
                paste0("Sensitivity = ", .)


        # plot as dot plot, x axis is cohort, y axis is y_pred, color is y_test, facet by ichorcna_strat
        data <- ans2 %>%
                filter(ichorcna_strat == target_ichorcna_strat)

        # add a col called cohort_detected_by_99threshold, if y_pred > target_99threshold, then cohort_detected_by_99threshold is 1, else 0
        data <- data %>%
                mutate(cohort_detected_by_99threshold = ifelse(y_pred > target_99threshold, 1, 0))
        # add a col called cohort_detected_by_95threshold, if y_pred > target_99threshold, then cohort_detected_by_99threshold is 1, else 0
        data <- data %>%
                mutate(cohort_detected_by_95threshold = ifelse(y_pred > target_95threshold, 1, 0))

        return(data)
}

# run the function recursively using ichorcna_strat_param_vec, using reduce
result <- lapply(ichorcna_strat_param_vec, FUN = annotate_detected, thresholds_param = thresholds, ans2_param = ans2, cnn_model_chosen_param = cnn_model_chosen)

# bind the result
result_bind <- bind_rows(result)


# save result_bind to a csv
result_bind %>%
        write_csv(cnn_raw_csv)


cnn_raw <- result_bind

###############################################################################
# parse the xgb data
###############################################################################

# get indtest scores

fl_indtest <- list.files(
        path = wd_xgb,
        "repeat_0_fold_0_metrics.csv", recursive = TRUE, full.names = TRUE
)

fl_indtest_read <- bettermc::mclapply(fl_indtest, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_indtest_read) <- fl_indtest

fl_indtest_read_bind <- bind_rows(fl_indtest_read, .id = "id") %>%
        mutate(ichorcna_strat = case_when(
                stringr::str_detect(id, "0-0\\.03") ~ "[0, 0.03]",
                stringr::str_detect(id, "0\\.03-0\\.1") ~ "(0.03, 0.1]",
                stringr::str_detect(id, "0\\.1-1") ~ "(0.1, 1]",
                stringr::str_detect(id, "all") ~ "all",
                TRUE ~ "unknown"
        ))
# set order of ichorcna_strat
fl_indtest_read_bind$ichorcna_strat <- factor(fl_indtest_read_bind$ichorcna_strat,
        levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all")
)

# make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
fl_indtest_read_bind <- fl_indtest_read_bind %>%
        mutate(feat = stringr::str_extract(id, "(?<=\\.feat\\.).*(?=\\.tf\\.)"))

# filter out feat contains xgb_all_feat_name
fl_indtest_read_bind <- fl_indtest_read_bind %>%
        filter(str_detect(feat, xgb_all_feat_name))


# plot the feature performance of 'cnv_sd_length_motif'
fl_indtest_read_bind <- fl_indtest_read_bind %>%
        # test_sen_98spe is 1 when test_auroc is 1
        # mutate(test_sen_98spe = ifelse(test_auroc == 1, 1, test_sen_98spe)) %>%
        # pivot_longer to make .metric as a column, all cols start with 'test_' to ".metric", values to ".estimate"
        pivot_longer(cols = starts_with("test_"), names_to = ".metric", values_to = ".estimate")

fl_indtest_read_bind_all_feature_model <- fl_indtest_read_bind %>%
        # filter(feat == layer_combination_chosen) %>%
        filter(model_name == xgb_model_chosen)

test_99spec_thresholds <- fl_indtest_read_bind_all_feature_model %>%
        filter(.metric == "test_threshold_99spe") %>%
        select(ichorcna_strat, .metric, .estimate) %>%
        distinct()

test_95spec_thresholds <- fl_indtest_read_bind_all_feature_model %>%
        filter(.metric == "test_threshold_95spe") %>%
        select(ichorcna_strat, .metric, .estimate) %>%
        distinct()

thresholds <- fl_indtest_read_bind_all_feature_model %>%
        select(ichorcna_strat, .metric, .estimate) %>%
        distinct()

###############################################################################

# meta_data <- read_csv(meta_csv_latest)

fl <- list.files(
        path = wd_xgb,
        "repeat_0_fold_0_y_test_y_pred_z_test.csv", recursive = TRUE, full.names = TRUE
)

fl_read <- bettermc::mclapply(fl, read_csv, show_col_types = FALSE, mc.cores = parallel::detectCores() / 2)

names(fl_read) <- fl

fl_read_bind <- bind_rows(fl_read, .id = "id") %>%
        mutate(ichorcna_strat = case_when(
                stringr::str_detect(id, "0-0\\.03") ~ "[0, 0.03]",
                stringr::str_detect(id, "0\\.03-0\\.1") ~ "(0.03, 0.1]",
                stringr::str_detect(id, "0\\.1-1") ~ "(0.1, 1]",
                stringr::str_detect(id, "all") ~ "all",
                TRUE ~ "unknown"
        ))
# set order of ichorcna_strat
fl_read_bind$ichorcna_strat <- factor(fl_read_bind$ichorcna_strat,
        levels = c("[0, 0.03]", "(0.03, 0.1]", "(0.1, 1]", "all")
)

# make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
fl_read_bind <- fl_read_bind %>%
        mutate(feat = stringr::str_extract(id, "(?<=\\.feat\\.).*(?=\\.tf\\.)"))

fl_read_bind <- fl_read_bind %>%
        filter(str_detect(feat, xgb_all_feat_name))

# # make a column called "feat" based on the path of the file, the feat should be everything between ".feat." and ".all_roc_curve_csv"
# fl_indtest_read_bind <- fl_indtest_read_bind %>%
#         mutate(feat = stringr::str_extract(id, "(?<=\\.feat\\.).*(?=\\.tf\\.)"))

# # filter out feat contains xgb_all_feat_name
# fl_indtest_read_bind <- fl_indtest_read_bind %>%
#         filter(str_detect(feat, xgb_all_feat_name))


# save fl_read_bind to a csv
# fl_read_bind %>%
#         write_csv(file.path(dirname(plot_file), "cnn_fl_scores_read_bind_indtest.csv"))


# only keep the all_feature model, by str_detect(id, "cnv_sd_length_ctRatio_slRatio")
fl_read_bind <- fl_read_bind %>%
        # filter(str_detect(feat, layer_combination_chosen)) %>%
        filter(model_name == xgb_model_chosen)


# # add meta_data$cohort to ans2, based on the bam_id
# meta_data_selected <- meta_data %>%
#         select(bam_id, cohort)

ans2 <- left_join(fl_read_bind, meta_data_selected, by = c("bam_id_test" = "bam_id"))

# y_test as factor
ans2$y_test <- factor(ans2$y_test, levels = c(0, 1))

# recode the y_test to "Healthy" and "Cancer"
ans2$y_test <- recode(ans2$y_test, "0" = "Healthy", "1" = "Cancer")

# make the cohort as factor
ans2$cohort <- factor(ans2$cohort)

# set the level of cohort, make "Healthy" the first level using relevel
ans2$cohort <- relevel(ans2$cohort, ref = "Healthy")

# report the number of samples in each ichorcna_strat and cohort
ans2_summary_xgb <- ans2 %>%
        group_by(ichorcna_strat, cohort) %>%
        summarize(n = n())

# stop running if ans2_summary_xgb and ans2_summary_cnn are not equal 

if (!isTRUE(all.equal(ans2_summary_xgb, ans2_summary_cnn))) {
  stop("The tibbles ans2_summary_xgb and ans2_summary_cnn are not the same.")
}

ans2 %>% group_by(ichorcna_strat, y_test) %>% summarize(n = n())

df <- expand.grid(ichorcna_strat_param = unique(ans2$ichorcna_strat))


ichorcna_strat_param_vec <- df$ichorcna_strat_param
# ans2


annotate_detected <- function(ichorcna_strat_param,
                              thresholds_param,
                              ans2_param,
                              cnn_model_chosen_param) {
        # plot
        target_ichorcna_strat <- ichorcna_strat_param
        target_99threshold <- thresholds %>%
                filter(.metric == "test_threshold_99spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)

        target_95threshold <- thresholds %>%
                filter(.metric == "test_threshold_95spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)

        sen_99spec <- thresholds %>%
                filter(.metric == "test_sen_99spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)

        sen_99_spec_label <- sen_99spec %>%
                round(., 4) %>%
                # convert to percentage and append %
                scales::percent(., accuracy = 0.01) %>%
                paste0("Sensitivity = ", .)
        sen_95spec <- thresholds %>%
                filter(.metric == "test_sen_95spe") %>%
                filter(ichorcna_strat == target_ichorcna_strat) %>%
                pull(.estimate)
        sen_95_spec_label <- sen_95spec %>%
                round(., 4) %>%
                # convert to percentage and append %
                scales::percent(., accuracy = 0.01) %>%
                paste0("Sensitivity = ", .)


        # plot as dot plot, x axis is cohort, y axis is y_pred, color is y_test, facet by ichorcna_strat
        data <- ans2 %>%
                filter(ichorcna_strat == target_ichorcna_strat)

        # add a col called cohort_detected_by_99threshold, if y_pred > target_99threshold, then cohort_detected_by_99threshold is 1, else 0
        data <- data %>%
                mutate(cohort_detected_by_99threshold = ifelse(y_pred > target_99threshold, 1, 0))
        # add a col called cohort_detected_by_95threshold, if y_pred > target_99threshold, then cohort_detected_by_99threshold is 1, else 0
        data <- data %>%
                mutate(cohort_detected_by_95threshold = ifelse(y_pred > target_95threshold, 1, 0))

        return(data)
}

# run the function recursively using ichorcna_strat_param_vec, using reduce
result <- lapply(ichorcna_strat_param_vec, FUN = annotate_detected, thresholds_param = thresholds, ans2_param = ans2, cnn_model_chosen_param = xgb_model_chosen)

# bind the result
result_bind <- bind_rows(result)


# save result_bind to a csv
result_bind %>%
        write_csv(xgb_raw_csv)


xgb_raw <- result_bind


###############################################################################
# plot the venn diagram
###############################################################################
xgb_cnn <- bind_rows(
        xgb_raw,
        cnn_raw
) %>%
        group_by(ichorcna_strat)

xgb_cnn %>%
        group_by(ichorcna_strat, model_name) %>%
        summarise(n = n())

xgb_cnn_split <- xgb_cnn %>%
        group_by(ichorcna_strat) %>%
        group_split(ichorcna_strat)


plot_venn <- function(x) {
        ichorcna_strat <- unique(x$ichorcna_strat)
        ichorcna_strat_label <- ichorcna_strat %>%
                # replace all blank with _
                stringr::str_replace_all(., " ", "_")
        # subsititute "resnet" with "ResNet" in cnn_model_chosen
        cnn_model_chosen_label <- cnn_model_chosen %>%
                stringr::str_replace_all(., "resnet", "ResNet")


        all_c <- x %>%
                filter(y_test == "Cancer") %>%
                pull(bam_id_test) %>%
                unique()
        all_h <- x %>%
                filter(y_test == "Healthy") %>%
                pull(bam_id_test) %>%
                unique()

        xgb_detected_95spe <- x %>%
                filter(model_name == "xgboost") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_95threshold == 1) %>%
                pull(bam_id_test) %>%
                unique()

        xgb_neg_detected_95spe <- x %>%
                filter(model_name == "xgboost") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_95threshold == 0) %>%
                pull(bam_id_test) %>%
                unique()

        xgb_detected_99spe <- x %>%
                filter(model_name == "xgboost") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_99threshold == 1) %>%
                pull(bam_id_test) %>%
                unique()

        xgb_neg_detected_99spe <- x %>%
                filter(model_name == "xgboost") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_99threshold == 0) %>%
                pull(bam_id_test) %>%
                unique()

        cnn_detected_95spe <- x %>%
                filter(model_name == "resnet18") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_95threshold == 1) %>%
                pull(bam_id_test) %>%
                unique()
        cnn_neg_detected_95spe <- x %>%
                filter(model_name == "resnet18") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_95threshold == 0) %>%
                pull(bam_id_test) %>%
                unique()
        cnn_detected_99spe <- x %>%
                filter(model_name == "resnet18") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_99threshold == 1) %>%
                pull(bam_id_test) %>%
                unique()
        cnn_neg_detected_99spe <- x %>%
                filter(model_name == "resnet18") %>%
                # filter(y_test == "Cancer") %>%
                filter(cohort_detected_by_99threshold == 0) %>%
                pull(bam_id_test) %>%
                unique()

        # venn diagram using the package ggVennDiagram
        ans_95 <- list(
                Cancer = all_c,
                xgb_95spe = xgb_detected_95spe,
                cnn_95spe = cnn_detected_95spe
        )

        # convert ans_95 to a dataframe
        # ans_95_df <- tibble(
        #         Cancer = all_c,
        #         xgb_95spe = xgb_detected_95spe,
        #         cnn_95spe = cnn_detected_95spe
        # )

        # # save ans_95_df to a csv
        # ans_95_df %>%
        #         write_csv(file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_95spe.csv")))

        ans_neg_95 <- list(
                Healthy = all_h,
                xgb_95spe = xgb_neg_detected_95spe,
                cnn_95spe = cnn_neg_detected_95spe
        )

        # ans_neg_95_df <- tibble(
        #         Healthy = all_h,
        #         xgb_95spe = xgb_neg_detected_95spe,
        #         cnn_95spe = cnn_neg_detected_95spe
        # )

        # # save ans_neg_95_df to a csv
        # ans_neg_95_df %>%
        #         write_csv(file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_neg_95spe.csv")))


        ans_99 <- list(
                Cancer = all_c,
                xgb_99spe = xgb_detected_99spe,
                cnn_99spe = cnn_detected_99spe
        )

        # ans_99_df <- tibble(
        #         Cancer = all_c,
        #         xgb_99spe = xgb_detected_99spe,
        #         cnn_99spe = cnn_detected_99spe
        # )

        # # save ans_99_df to a csv
        # ans_99_df %>%
        #         write_csv(file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_99spe.csv")))

        ans_neg_99 <- list(
                Healthy = all_h,
                xgb_99spe = xgb_neg_detected_99spe,
                cnn_99spe = cnn_neg_detected_99spe
        )

        # ans_neg_99_df <- tibble(
        #         Healthy = all_h,
        #         xgb_99spe = xgb_neg_detected_99spe,
        #         cnn_99spe = cnn_neg_detected_99spe
        # )

        # # save ans_neg_99_df to a csv
        # ans_neg_99_df %>%
        #         write_csv(file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_neg_99spe.csv")))


        # ans_99 --------------------------------------------------------------
        p99 <- ggVennDiagram(ans_99,
                category.names = c("True Pos", "XGBoost", cnn_model_chosen_label),
                label_alpha = 0,
                label = "count",
                set_color = c("black", "darkblue", "orange3"),
                set_size = c(3, 3, 3),
                label_size = 3,
                edge_size = 1
        )

        p99 <- p99 +
                scale_x_continuous(expand = expansion(mult = .2)) +
                labs(
                        title = "   ",
                        subtitle = "   "
                ) +
                # title text size to 6
                theme(plot.title = element_text(size = 7)) +
                # subtitle text size to 6
                theme(plot.subtitle = element_text(size = 7)) +
                # remove the legend
                theme(legend.position = "none") +
                scale_fill_gradient(low = "white", high = "grey90")


        ggsave(
                filename = file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_99spe.pdf")),
                plot = p99,
                width = 6,
                height = 6,
                units = "cm"
        )
        message("Saved to: ", file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_99spe.pdf")))

        # ans_neg_99 --------------------------------------------------------------
        p_neg99 <- ggVennDiagram(ans_neg_99,
                category.names = c("True Neg", "XGBoost", cnn_model_chosen_label),
                label_alpha = 0,
                label = "count",
                set_color = c("black", "darkblue", "orange3"),
                set_size = c(3, 3, 3),
                label_size = 3,
                edge_size = 1
        )

        p_neg99 <- p_neg99 +
                scale_x_continuous(expand = expansion(mult = .2)) +
                labs(
                        title = "   ",
                        subtitle = "   "
                ) +
                # title text size to 6
                theme(plot.title = element_text(size = 7)) +
                # subtitle text size to 6
                theme(plot.subtitle = element_text(size = 7)) +
                # remove the legend
                theme(legend.position = "none") +
                scale_fill_gradient(low = "white", high = "grey90")

        ichorcna_strat_label <- ichorcna_strat %>%
                # replace all blank with _
                stringr::str_replace_all(., " ", "_")

        ggsave(
                filename = file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_neg_99spe.pdf")),
                plot = p_neg99,
                width = 6,
                height = 6,
                units = "cm"
        )
        message("Saved to: ", file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_neg_99spe.pdf")))




        # ans_95 --------------------------------------------------------------
        p95 <- ggVennDiagram(ans_95,
                category.names = c("True Pos", "XGBoost", cnn_model_chosen_label),
                label_alpha = 0,
                label = "count",
                set_color = c("black", "darkblue", "orange3"),
                set_size = c(3, 3, 3),
                label_size = 3,
                edge_size = 1
        )

        p95 <- p95 +
                scale_x_continuous(expand = expansion(mult = .2)) +
                labs(
                        title = "   ",
                        subtitle = "   "
                ) +
                # title text size to 6
                theme(plot.title = element_text(size = 7)) +
                # subtitle text size to 6
                theme(plot.subtitle = element_text(size = 7)) +
                # remove the legend
                theme(legend.position = "none") +
                scale_fill_gradient(low = "white", high = "grey90")


        ichorcna_strat_label <- ichorcna_strat %>%
                # replace all blank with _
                stringr::str_replace_all(., " ", "_")

        ggsave(
                filename = file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_95spe.pdf")),
                plot = p95,
                width = 6,
                height = 6,
                units = "cm"
        )
        message("Saved to: ", file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_95spe.pdf")))

        # ans_neg_95 --------------------------------------------------------------
        p_neg95 <- ggVennDiagram(ans_neg_95,
                category.names = c("True Neg", "XGBoost", cnn_model_chosen_label),
                label_alpha = 0,
                label = "count",
                set_color = c("black", "darkblue", "orange3"),
                set_size = c(3, 3, 3),
                label_size = 3,
                edge_size = 1
        )

        p_neg95 <- p_neg95 +
                scale_x_continuous(expand = expansion(mult = .2)) +
                labs(
                        title = "   ",
                        subtitle = "   "
                ) +
                # title text size to 6
                theme(plot.title = element_text(size = 7)) +
                # subtitle text size to 6
                theme(plot.subtitle = element_text(size = 7)) +
                # remove the legend
                theme(legend.position = "none") +
                scale_fill_gradient(low = "white", high = "grey90")

        ichorcna_strat_label <- ichorcna_strat %>%
                # replace all blank with _
                stringr::str_replace_all(., " ", "_")

        ggsave(
                filename = file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_neg_95spe.pdf")),
                plot = p_neg95,
                width = 6,
                height = 6,
                units = "cm"
        )
        message("Saved to: ", file.path(plot_path, paste0("venn_", ichorcna_strat_label, "_neg_95spe.pdf")))
}

lapply(xgb_cnn_split, plot_venn)



# save all to venn.rda
save.image(file = file.path(plot_path, "venn_analysis.rda"))



###############################################################################
# plot the tumor fraction of cancer samples
###############################################################################

# get the cols needed from
meta <- read_csv(meta_csv_latest)
meta_clean <- meta %>%
        select(bam_id, clinical_tf, ichorcna_tf)

xgb_cnn_annotated <- left_join(xgb_cnn, meta_clean, by = c("bam_id_test" = "bam_id"))

xgb_cnn_annotated_filtered <- xgb_cnn_annotated %>%
        filter(y_test == "Cancer") %>%
        # pivot longer, gather "cohort_detected_by_99threshold", "cohort_detected_by_95threshold"
        pivot_longer(
                cols = c("cohort_detected_by_99threshold", "cohort_detected_by_95threshold"),
                names_to = "detect_by", values_to = "detected"
        )




xgb_95 <- xgb_cnn_annotated_filtered %>%
        filter(detect_by == "cohort_detected_by_95threshold") %>%
        filter(model_name == "xgboost")

cnn_95 <- xgb_cnn_annotated_filtered %>%
        filter(detect_by == "cohort_detected_by_95threshold") %>%
        filter(model_name == "resnet18")

xgb_99 <- xgb_cnn_annotated_filtered %>%
        filter(detect_by == "cohort_detected_by_99threshold") %>%
        filter(model_name == "xgboost")

cnn_99 <- xgb_cnn_annotated_filtered %>%
        filter(detect_by == "cohort_detected_by_99threshold") %>%
        filter(model_name == "resnet18")


run_list <- list(
        xgb_95,
        cnn_95,
        xgb_99,
        cnn_99
)

tf_plot_func <- function(ichorcna_strat_param, x) {
        xraw <- x
        x <- xraw %>%
                filter(ichorcna_strat == ichorcna_strat_param) %>%
                group_by(detected) %>%
                arrange(desc(ichorcna_tf)) %>%
                ungroup()

        x$bam_id_test <- factor(x$bam_id_test, levels = x$bam_id_test)

        model_name <- x %>%
                pull(model_name) %>%
                unique()

        detect_by <- x %>%
                pull(detect_by) %>%
                unique()


        ichorcna_strat_param_label <- ichorcna_strat_param %>%
                # replace all blank with _
                stringr::str_replace_all(., " ", "_")

        tf_plotfile <- file.path(
                plot_path,
                paste("ichore_tumor_fraction_plot", ichorcna_strat_param_label, model_name, detect_by, ".pdf", sep = "_")
        )

        # set x$detected as factor
        x <- x %>%
                mutate(detected_label = ifelse(detected == 1, "Detected", "Not Detected")) %>%
                group_by(detected_label) %>%
                mutate(detected_count = n()) %>%
                mutate(detected_facet_label = paste0(detected_label, " (n = ", detected_count, ")")) %>%
                ungroup()

        detected_facet_label_vec <- x %>%
                pull(detected_facet_label) %>%
                unique()

        detected_facet_label_vec_color <- c("red", "lightgrey")
        names(detected_facet_label_vec_color) <- detected_facet_label_vec




        x$detected <- factor(x$detected)

        n_cancer_samples <- x %>%
                select(bam_id_test) %>%
                unique() %>%
                nrow()
        x_lab <- paste0("Cancer Samples (n = ", n_cancer_samples, ")")


        p_tf <- ggplot(x, aes(x = bam_id_test, y = ichorcna_tf)) +
                geom_point(aes(x = bam_id_test, y = ichorcna_tf, fill = detected_facet_label),
                        size = 0.5,
                        shape = 21,
                        stroke = 0.1
                ) +
                labs(
                        x = x_lab,
                        y = "ichorCNA Tumor Fraction"
                ) +
                theme_classic() +
                # set color manually
                scale_color_manual(values = detected_facet_label_vec_color) +
                theme(
                        legend.position = "none",
                        axis.text.x = element_blank()
                ) +
                facet_wrap(~detected_facet_label, scales = "free_x") +
                # make strip box height smaller
                theme(strip.text = element_text(size = 6)) +
                # remove x axis ticks
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
                # remove x axis label
                theme(axis.title.x = element_blank()) +
                # set strip panel stroke color to none
                theme(strip.background = element_rect(color = NA)) +
                # set strip panel color to lightgrey
                theme(strip.background = element_rect(fill = "lightgrey")) +
                # y title size to 6
                theme(axis.title.y = element_text(size = 6))
        ggsave(
                filename = tf_plotfile,
                plot = p_tf,
                width = 9,
                height = 3.5,
                units = "cm"
        )
        message("Saved to: ", tf_plotfile)


        #--------------------------------------------------------------------------
        # clinical tf plot
        #--------------------------------------------------------------------------

        x <- xraw %>%
                filter(ichorcna_strat == ichorcna_strat_param) %>%
                group_by(detected) %>%
                arrange(desc(clinical_tf)) %>%
                ungroup()

        x$bam_id_test <- factor(x$bam_id_test, levels = x$bam_id_test)

        model_name <- x %>%
                pull(model_name) %>%
                unique()

        detect_by <- x %>%
                pull(detect_by) %>%
                unique()


        ichorcna_strat_param_label <- ichorcna_strat_param %>%
                # replace all blank with _
                stringr::str_replace_all(., " ", "_")

        tf_plotfile <- file.path(
                plot_path,
                paste("clinical_tumor_fraction_plot", ichorcna_strat_param_label, model_name, detect_by, ".pdf", sep = "_")
        )

        # set x$detected as factor
        x <- x %>%
                mutate(detected_label = ifelse(detected == 1, "Detected", "Not Detected")) %>%
                group_by(detected_label) %>%
                mutate(detected_count = n()) %>%
                mutate(detected_facet_label = paste0(detected_label, " (n = ", detected_count, ")")) %>%
                ungroup()

        detected_facet_label_vec <- x %>%
                pull(detected_facet_label) %>%
                unique()

        detected_facet_label_vec_color <- c("#dd1d1d", "lightgrey")
        names(detected_facet_label_vec_color) <- detected_facet_label_vec




        x$detected <- factor(x$detected)

        n_cancer_samples <- x %>%
                select(bam_id_test) %>%
                unique() %>%
                nrow()
        x_lab <- paste0("Cancer Samples (n = ", n_cancer_samples, ")")

        # make x2 , chaging the clinical_tf to clinical_tf2, NA to "NA"
        x$clinical_tf2 <- ifelse(is.na(x$clinical_tf), -0.1, x$clinical_tf)

        # get the min clinical_tf in the detected == 1
        min_clinical_tf <- x %>%
                filter(detected == 1) %>%
                # filter out NA
                filter(!is.na(clinical_tf)) %>%
                pull(clinical_tf) %>%
                min()

        print(x %>%
                filter(detected == 1) %>%
                # filter out NA
                filter(!is.na(clinical_tf)) %>%
                pull(clinical_tf))

        anno_text <- paste0("LOD: ", min_clinical_tf)


        p_tf <- ggplot(x, aes(x = bam_id_test, y = clinical_tf2)) +
                geom_point(aes(x = bam_id_test, y = clinical_tf2, fill = detected_facet_label),
                        size = 1,
                        shape = 21,
                        stroke = 0.1
                ) +
                # geom_point(aes(y = clinical_tf2),
                #         color = "black",
                #         data = x %>% filter(is.na(clinical_tf)),
                #         size = 0.3,
                #         shape = 21,
                #         fill = "white",
                #         stroke = 0.1
                # ) +
                # add text annotation
                annotate("text", x = 1, y = min_clinical_tf, label = anno_text, color = "black", size = 2) +
                labs(
                        x = x_lab,
                        y = "Clinically Validated Tumor Fraction"
                ) +
                theme_classic() +
                # change the y axis label to "NA" when it is -1
                scale_y_continuous(labels = function(x) ifelse(x == -0.1, "NA", x)) +
                # set color manually
                scale_color_manual(values = detected_facet_label_vec_color) +
                theme(
                        legend.position = "none",
                        axis.text.x = element_blank()
                ) +
                facet_wrap(~detected_facet_label, scales = "free_x") +
                # make strip box height smaller
                theme(strip.text = element_text(size = 6)) +
                # remove x axis ticks
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
                # remove x axis label
                theme(axis.title.x = element_blank()) +
                # set strip panel stroke color to none
                theme(strip.background = element_rect(color = NA)) +
                # set strip panel color to lightgrey
                theme(strip.background = element_rect(fill = "lightgrey")) +
                # y title size to 6
                theme(axis.title.y = element_text(size = 5))
        ggsave(
                filename = tf_plotfile,
                plot = p_tf,
                width = 9,
                height = 3.5,
                units = "cm"
        )
        message("Saved to: ", tf_plotfile)
}

for (x in run_list) {
        lapply(ichorcna_strat_param_vec, FUN = tf_plot_func, x = x)
}
