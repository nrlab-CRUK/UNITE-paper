# library(ggupset)
# library(ggplot2)
# library(tidyverse, warn.conflicts = FALSE)
# library(argparse)
# library(kableExtra)
# library(gridExtra)
# library(grid)
# library(gghalves)
# library(nord)
# library(ggrepel)
# library(gghalves) and install it if not exists
# if (!requireNamespace("gghalves", quietly = TRUE)) {
#         install.packages("gghalves")
# }

if (!require("pacman")) install.packages("pacman")

# Load multiple packages
pacman::p_load(tidyverse,
ggupset,
gghalves,
ggrepel,
ggplot2,
argparse,
kableExtra,
gridExtra,
grid,
patchwork, 
MultiAssayExperiment,
nord,
plotly,
ggmagnify,
gghighlight,
ggdist,
htmlwidgets,
ggpubr)


source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

# argparse
parser <- argparse::ArgumentParser()
parser$add_argument("--cnn_data_path", type = "character", help = "cnn_data_path", default = "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn14/")
parser$add_argument("--xgb_data_path", type = "character", help = "xgb_data_path", default = "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb10/")
parser$add_argument("--cnn_layers", type = "character", help = "cnn_layers", default = "C4")
parser$add_argument("--xgb_feat", type = "character", help = "xgb_feat", default = "feat.cnv_sd_length_ctRatio_slRatio")
parser$add_argument("--cnn_model_name", type = "character", help = "cnn_model_name", default = "resnet18")
parser$add_argument("--tf_min", type = "character", help = "tf_min", default = "0")
parser$add_argument("--tf_max", type = "character", help = "tf_max", default = "1")
args <- parser$parse_args()
cnn_data_path <- args$cnn_data_path
xgb_data_path <- args$xgb_data_path
cnn_layers <- args$cnn_layers
xgb_feat <- args$xgb_feat
cnn_model_name <- args$cnn_model_name
tf_min <- args$tf_min
tf_max <- args$tf_max
xgb_tf_label <- paste0(tf_min, "-", tf_max, "_pairwise")
cnn_tf_label <- paste0("tf_", tf_min, "_", tf_max)

# if tf_min ==0 and tf_max ==1, set xgb_tf_label to "all"
if (tf_min == "0" && tf_max == "1") {
        xgb_tf_label <- "all"
}


xgb_metric_file <- list.files(xgb_data_path, pattern = "repeat_0_fold_0_metrics.csv", recursive = TRUE, full.names = TRUE) %>%
        as_tibble() %>%
        filter(str_detect(value, xgb_tf_label)) %>%
        filter(str_detect(value, xgb_feat)) %>%
        pull()

# stop if xgb_metric_file is empty or contains more than 1 element
if (length(xgb_metric_file) != 1) {
        stop("xgb_metric_file is empty or contains more than 1 element")
}

cnn_metric_file <- list.files(cnn_data_path, pattern = "repeat_0_fold_0_metrics.csv", recursive = TRUE, full.names = TRUE) %>%
        as_tibble() %>%
        filter(str_detect(value, cnn_layers)) %>%
        filter(str_detect(value, cnn_model_name)) %>%
        filter(str_detect(value, cnn_tf_label)) %>%
        # filter out values containing "mismatch"
        filter(!str_detect(value, "mismatch")) %>%
        pull()

# stop if cnn_metric_file is empty or contains more than 1 element
if (length(cnn_metric_file) != 1) {
        stop("cnn_metric_file is empty or contains more than 1 element")
}


xgboost_roc_file <- list.files(xgb_data_path, pattern = "repeat_0_fold_0_roc.csv", recursive = TRUE, full.names = TRUE) %>%
        as_tibble() %>%
        filter(str_detect(value, xgb_tf_label)) %>%
        filter(str_detect(value, xgb_feat)) %>%
        pull()

# stop if xgboost_roc_file is empty or contains more than 1 element
if (length(xgboost_roc_file) != 1) {
        stop("xgboost_roc_file is empty or contains more than 1 element")
}

cnn_roc_file <- list.files(cnn_data_path, pattern = "repeat_0_fold_0_roc.csv", recursive = TRUE, full.names = TRUE) %>%
        as_tibble() %>%
        filter(str_detect(value, cnn_layers)) %>%
        filter(str_detect(value, cnn_model_name)) %>%
        filter(str_detect(value, cnn_tf_label)) %>%
        # filter out values containing "mismatch"
        filter(!str_detect(value, "mismatch")) %>%
        pull()

# stop if cnn_roc_file is empty or contains more than 1 element
if (length(cnn_roc_file) != 1) {
        stop("cnn_roc_file is empty or contains more than 1 element")
}
# xgb_metric_file <- "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb9/0-0.03_pairwise/xgboost_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0-0.03_pairwise_roc_curve_csv/ind_test/repeat_0_fold_0_metrics.csv"
# cnn_metric_file <- "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn14/plasma_tp1_0.1x_model_C4/roc_curve_resnet18_tf_0_0.03/ind_test_random/resnet18/repeat_0_fold_0_metrics.csv"

# xgboost_roc_file <- "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb9/0-0.03_pairwise/xgboost_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0-0.03_pairwise_roc_curve_csv/ind_test/repeat_0_fold_0_roc.csv"
# cnn_roc_file <- "/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn14/plasma_tp1_0.1x_model_C4/roc_curve_resnet18_tf_0_0.03/ind_test_random/resnet18/repeat_0_fold_0_roc.csv"
plot_file <- file.path(
        "/home/nrlab/wang04/ulyses/4_ROC/", cnn_layers,
        paste0("tf_", tf_min, "_", tf_max, "_ind_test_roc.pdf")
)



best_xgboost_roc <- read_csv(xgboost_roc_file) %>%
        mutate(model = "XGBoost")


best_cnn_roc <- read_csv(cnn_roc_file) %>%
        mutate(model = "ResNet18")

xgb_auroc <- read_csv(xgb_metric_file) %>%
        pull("test_auroc")

xgb_threshold99spec <- read_csv(xgb_metric_file) %>%
        pull("test_threshold_99spe")

xgb_threshold95spec <- read_csv(xgb_metric_file) %>%
        pull("test_threshold_95spe")

cnn_auroc <- read_csv(cnn_metric_file) %>%
        pull("test_auroc")

cnn_threshold99spec <- read_csv(cnn_metric_file) %>%
        pull("test_threshold_99spe")

cnn_threshold95spec <- read_csv(cnn_metric_file) %>%
        pull("test_threshold_95spe")

all <- rbind(best_xgboost_roc, best_cnn_roc)

all <- all %>%
        mutate(model_label = ifelse(model == "XGBoost", paste("XGBoost (AUROC = ", round(xgb_auroc, 3), ")", sep = ""), paste("ResNet18 (AUROC = ", round(cnn_auroc, 3), ")", sep = "")))

m1 <- all$model_label |>
        unique() |>
        pluck(1)
m2 <- all$model_label |>
        unique() |>
        pluck(2)

all$model_label <- factor(all$model_label, levels = c(m1, m2))
colors <- c("darkblue", "orange3") %>%
        set_names(c(m1, m2))




# plot roc curve, group by model
roc_plot <- ggplot(all, aes(x = fpr, y = tpr, color = model_label, group = model)) +
        geom_line(linewidth = 0.5) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#bbb9b9", linewidth = 0.1) +
        theme_classic() +
        theme_NPJ() +
        # set legend position
        theme(legend.position = "inside") +
        scale_color_manual(values = colors) +
        labs(x = "1 - Specificity", y = "Sensitivity") +
        theme(legend.position.inside = c(0.7, 0.3)) +
        # reduce the distance between legend elements
        theme(legend.key.size = unit(0.2, "cm")) +
        theme(legend.title = element_blank()) +
        theme(legend.text = element_text(size = 5.5)) +
        theme(axis.text = element_text(size = 6)) +
        theme(axis.title = element_text(size = 7))

# # add a point to the roc_plot where the sensitivity is 0.99
# roc_plot <- roc_plot +
#         geom_point(data = all %>% filter(thresholds > xgb_threshold99spec) %>% filter(model == "XGBoost") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr), color = "darkblue", size = 2, shape = 21) +
#         geom_point(data = all %>% filter(thresholds > cnn_threshold99spec) %>% filter(model == "ResNet18") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr), color = "orange3", size = 2, shape = 21) +
#         geom_text_repel(data = all %>% filter(thresholds > xgb_threshold99spec) %>% filter(model == "XGBoost") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr, label = paste("99% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01))), hjust = -0.1, vjust = 1, size = 2.5, color = "darkblue") +
#         geom_text_repel(data = all %>% filter(thresholds > cnn_threshold99spec) %>% filter(model == "ResNet18") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr, label = paste("99% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01))), hjust = -0.1, vjust = 0, size = 2.5, color = "orange3")

# # add point to the roc_plot where the specificity is 0.95

# roc_plot <- roc_plot +
#         geom_point(data = all %>% filter(thresholds > xgb_threshold95spec) %>% filter(model == "XGBoost") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr), color = "darkblue", size = 2, shape = 21) +
#         geom_point(data = all %>% filter(thresholds > cnn_threshold95spec) %>% filter(model == "ResNet18") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr), color = "orange3", size = 2, shape = 21) +
#         geom_text_repel(data = all %>% filter(thresholds > xgb_threshold95spec) %>% filter(model == "XGBoost") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr, label = paste("95% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01))), hjust = -0.1, vjust = 1.5, size = 2.5, color = "darkblue") +
#         geom_text_repel(data = all %>% filter(thresholds > cnn_threshold95spec) %>% filter(model == "ResNet18") %>% slice_max(tpr, n = 1), aes(x = fpr, y = tpr, label = paste("95% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01))), hjust = -0.1, vjust = 0.3, size = 2.5, color = "orange3")


# simplify this
xgb_99 <- all %>%
        filter(thresholds > xgb_threshold99spec) %>%
        filter(model == "XGBoost") %>%
        slice_max(tpr, n = 1) %>%
        mutate(text = paste("99% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01)))

xgb_95 <- all %>%
        filter(thresholds > xgb_threshold95spec) %>%
        filter(model == "XGBoost") %>%
        slice_max(tpr, n = 1) %>%
        mutate(text = paste("95% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01)))
cnn_99 <- all %>%
        filter(thresholds > cnn_threshold99spec) %>%
        filter(model == "ResNet18") %>%
        slice_max(tpr, n = 1) %>%
        mutate(text = paste("99% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01)))

cnn_95 <- all %>%
        filter(thresholds > cnn_threshold95spec) %>%
        filter(model == "ResNet18") %>%
        slice_max(tpr, n = 1) %>%
        mutate(text = paste("95% Spe, Sen = ", scales::percent(tpr, accuracy = 0.01)))

annotation_tibble <- bind_rows(xgb_99, xgb_95, cnn_99, cnn_95)

roc_plot2 <- roc_plot +
        geom_point(
                data = annotation_tibble,
                aes(x = fpr, y = tpr, group = model_label, color = model_label),
                show.legend = FALSE,
                size = 2,
                shape = 21
        ) +
        geom_text_repel(
                data = annotation_tibble,
                aes(x = fpr, y = tpr, group = model_label, color = model_label, label = text),
                # color = "grey50",
                segment.color = "black",
                segment.size = 0.08,
                show.legend = FALSE,
                hjust = -0.1, vjust = 1.5, size = 2.5
        )


ggsave(plot_file, roc_plot2,
        width = 6, height = 6, units = "cm",
        dpi = 300,
        # create folder if not exists
        create.dir = TRUE
)

message("ROC curve plot saved to ", plot_file)
