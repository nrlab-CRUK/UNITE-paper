library(tidyverse)
library(gghalves)
library(ggpubr)
library(ggdist)

# convert which_tf_strat to a parameter
parser <- argparse::ArgumentParser()
parser$add_argument("--which_tf_strat", type = "character", help = "which_tf_strat", default = "all")
parser$add_argument("--cnn_model_to_keep", type = "character", help = "cnn_model_to_keep", default = "resnet18")
parser$add_argument("--xgb_model_to_keep", type = "character", help = "xgb_model_to_keep", default = "xgboost")
parser$add_argument("--lr_model_to_keep", type = "character", help = "lr_model_to_keep", default = "lr")
parser$add_argument("--feat_chosen", type = "character", help = "feat_chosen", default = c("C4", "All", "TF"))
parser$add_argument("--cnn_feat_chosen", type = "character", help = "feat_chosen", default = c("C4"))
parser$add_argument("--xgb_metrics_file", type = "character", help = "xgb_metrics_file", default = "/home/nrlab/wang04/ulyses/1_nested_cv_plots/xgboost_ranking_plot_main_fig.csv")
parser$add_argument("--cnn_metrics_file", type = "character", help = "cnn_metrics_file", default = "/home/nrlab/wang04/ulyses/1_nested_cv_plots/resnet18/C4/cnn_ranking_plot_main_fig.csv")
parser$add_argument("--lr_metrics_file", type = "character", help = "lr_metrics_file", default = "/home/nrlab/wang04/ulyses/1_nested_cv_plots/lr_ranking_plot_main_fig.csv")
# add a parameter called "which_metrics_to_plot"
parser$add_argument("--which_metrics_to_plot", type = "character", help = "which_metrics_to_plot", default = c(
        "AUROC",
        "AUPRC",
        "Sensitivity at 99% Specificity",
        "Sensitivity at 95% Specificity"
))
parser$add_argument("--all_hight_color", type = "character", help = "all_hight_color", default = "#c1dde7")
# parse the args
args <- parser$parse_args()
which_tf_strat <- args$which_tf_strat
cnn_model_to_keep <- args$cnn_model_to_keep
xgb_model_to_keep <- args$xgb_model_to_keep
lr_model_to_keep <- args$lr_model_to_keep
feat_chosen <- args$feat_chosen
cnn_feat_chosen <- args$cnn_feat_chosen
which_metrics_to_plot <- args$which_metrics_to_plot
all_hight_color <- args$all_hight_color

# metrics files and read in ----------------------------------------------------
xgb_metrics_file <- args$xgb_metrics_file
cnn_metrics_file <- args$cnn_metrics_file
lr_metrics_file <- args$lr_metrics_file
xgb_metrics <- read_csv(xgb_metrics_file)
cnn_metrics <- read_csv(cnn_metrics_file)
lr_metrics <- read_csv(lr_metrics_file)
# ------------------------------------------------------------------------------

#------------------------------plots parameters---------------------------------
manual_colors <- c("XGBoost" = "darkblue", "ResNet18" = "orange3", "LR" = "#2c2b2b")
pd <- position_dodge(0.4)
plot_unit <- "mm"
#-------------------------------------------------------------------------------

# preprocessing
source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
model_to_keep <- c(cnn_model_to_keep, xgb_model_to_keep, lr_model_to_keep)
y_axis_lab <- paste0("Performance at ", which_tf_strat)
plotfile <- file.path(
        "/home/nrlab/wang04/ulyses/1_nested_cv_plots/",
        "direct_comp_xgb_and_cnn",
        paste0("xgb_cnn_comparison_", which_tf_strat, ".pdf")
)


main_fig_cnn_vs_xgb_plot <- file.path(
        "/home/nrlab/wang04/ulyses/1_nested_cv_plots/",
        "direct_comp_xgb_and_cnn",
        paste0("xgb_vs_cnn_main_fig", ".pdf")
)

main_fig_between_cnns <- file.path(
        "/home/nrlab/wang04/ulyses/1_nested_cv_plots/",
        "direct_comp_xgb_and_cnn",
        paste0("cnn_models_internal_comparison", ".pdf")
)

# read in the data
xgb_cnn <- bind_rows(xgb_metrics, cnn_metrics, lr_metrics)
# set ichorcna_strat to factor, levels to ichorcna_tf_strat_hyper
xgb_cnn$ichorcna_strat <- factor(xgb_cnn$ichorcna_strat, levels = ichorcna_tf_strat_for_model_performance_hyper)


################################################################################
# TODO: main figure: xgb all VS cnn C5
################################################################################

# filter the data
xgb_cnn_filtered <- xgb_cnn %>%
        filter(.metric == "AUROC") %>%
        # filter(ichorcna_strat == which_tf_strat) %>%
        filter(model_name %in% model_to_keep) %>%
        filter(feat %in% feat_chosen) %>%
        mutate(model_name_label = case_when(
                model_name == "xgboost" ~ "XGBoost",
                model_name == "resnet18" ~ "ResNet18",
                model_name == "lr" ~ "LR",
                TRUE ~ model_name
        ))

# set the levels of model_name to the same as model_to_keep
xgb_cnn_filtered$model_name_label <- factor(xgb_cnn_filtered$model_name_label, levels = c("ResNet18", "XGBoost", "LR"))


p_main <- ggplot(data = xgb_cnn_filtered, aes(x = ichorcna_strat, y = .estimate, color = model_name_label, fill = model_name_label, group = model_name_label)) +
        # Add background color to 'all' group
        geom_rect(
                data = subset(xgb_cnn_filtered, ichorcna_strat == "all"),
                aes(
                        xmin = as.numeric(ichorcna_strat) - 0.5,
                        xmax = as.numeric(ichorcna_strat) + 0.5,
                        ymin = -Inf,
                        ymax = Inf
                ),
                fill = all_hight_color,
                alpha = 0.3,
                inherit.aes = FALSE
        ) +
        stat_compare_means(method = "kruskal.test", label = "p.signif", label.y = c(0.92, 0.991, 0.995, 0.930), size = 2.5, color = "#756a6a") +
        stat_summary(fun = mean, geom = "line", data = xgb_cnn_filtered %>% filter(ichorcna_strat != "all"), linewidth = 0.8) +
        stat_summary(fun = mean, geom = "point", shape = 21, color = "grey2", size = 2) +
        stat_summary(
                fun.data = mean_cl_boot,
                geom = "errorbar",
                position = pd,
                linewidth = 0.3,
                width = 0.25
        ) +
        facet_wrap(~.metric, scales = "free_y", drop = TRUE) +
        # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
        scale_color_manual(values = manual_colors) +
        scale_fill_manual(values = manual_colors) +
        theme_classic() +
        # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
        # add legend
        theme(strip.background = element_rect(fill = "#e9e6e6", color = "#e9e6e6")) +
        labs(x = NULL, y = "Mean performancce and 95% CI") +
        # remove x-axis coordiante extension
        scale_x_discrete(expand = c(0.05, 0)) +
        theme(panel.spacing = unit(0.1, "cm")) +
        theme(strip.text = element_text(size = 6)) +
        theme(axis.text = element_text(size = 6)) +
        # remove legend title
        theme(axis.title = element_text(size = 6)) +
        theme(legend.position = "right") +
        theme(legend.text = element_text(size = 6), legend.title = element_blank()) +
        theme(legend.key = element_rect(fill = "transparent", color = "transparent")) +
        theme(
                legend.box.background = element_rect(fill = "transparent", color = "transparent"),
                legend.box.spacing = unit(0, "pt"), # The spacing between the plotting area and the legend box (unit)
                legend.margin = margin(0, 0, 0, 0),
                plot.margin = margin(0, 0, 0, 0)
        )



ggsave(main_fig_cnn_vs_xgb_plot,
        p_main,
        width = 75,
        height = 60,
        unit = plot_unit,
        create.dir = TRUE
)

message(main_fig_cnn_vs_xgb_plot)

# save the plot to rds file
plot_file_rds <- gsub(".pdf", ".rds", main_fig_cnn_vs_xgb_plot)
saveRDS(p_main, plot_file_rds)

# hide legend
p_main_no_legend <- p_main + theme(legend.position = "none")
plot_file_no_legend <- gsub(".pdf", "_no_legend.pdf", main_fig_cnn_vs_xgb_plot)
ggsave(plot_file_no_legend,
        p_main_no_legend,
        width = 50,
        height = 60,
        unit = plot_unit,
        create.dir = TRUE
)
message(plot_file_no_legend)


################################################################################
# main figure: cnn models internal comparison
################################################################################

cnn_models_to_compare <- c(
        "resnet18",
        "resnet34",
        "tf_efficientnetv2_s.in1k",
        "plain_cnn",
        "vit_base"
)
manual_colors2 <- c(
        "ResNet18" = "orange3",
        "ResNet34" = "#d44bb2",
        "EfficientNetV2S" = "#566156",
        "Plain CNN" = "#2c892c",
        "ViT Base" = "#320162"
)
# filter the data
xgb_cnn_filtered <- xgb_cnn %>%
        filter(.metric == "AUROC") %>%
        filter(model_name %in% cnn_models_to_compare) %>%
        filter(feat %in% cnn_feat_chosen) %>%
        mutate(model_name_label = case_when(
                model_name == "resnet34" ~ "ResNet34",
                model_name == "resnet18" ~ "ResNet18",
                model_name == "tf_efficientnetv2_s.in1k" ~ "EfficientNetV2S",
                model_name == "plain_cnn" ~ "Plain CNN",
                model_name == "vit_base" ~ "ViT Base",
                TRUE ~ model_name
        ))

p_main <- ggplot(data = xgb_cnn_filtered, aes(x = ichorcna_strat, y = .estimate, color = model_name_label, fill = model_name_label, group = model_name_label)) +
        # Add background color to 'all' group
        geom_rect(
                data = subset(xgb_cnn_filtered, ichorcna_strat == "all"),
                aes(
                        xmin = as.numeric(ichorcna_strat) - 0.5,
                        xmax = as.numeric(ichorcna_strat) + 0.5,
                        ymin = -Inf,
                        ymax = Inf
                ),
                fill = all_hight_color,
                alpha = 0.3,
                inherit.aes = FALSE
        ) +
        stat_compare_means(
                method = "kruskal", label = "p.signif", label.y = c(0.91, 0.96, 1, 0.93), size = 3, color = "#756a6a",
                paired = FALSE
        ) +
        stat_summary(fun = mean, geom = "line", data = xgb_cnn_filtered %>% filter(ichorcna_strat != "all"), linewidth = 0.8) +
        stat_summary(fun = mean, geom = "point", shape = 21, color = "grey2", size = 2) +
        stat_summary(
                fun.data = mean_cl_boot,
                geom = "errorbar",
                position = pd,
                linewidth = 0.3,
                width = 0.25
        ) +
        facet_wrap(~.metric, scales = "free_y", drop = TRUE) +
        # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
        scale_color_manual(values = manual_colors2) +
        scale_fill_manual(values = manual_colors2) +
        theme_classic() +
        # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
        # add legend
        theme(strip.background = element_rect(fill = "#e9e6e6", color = "#e9e6e6")) +
        labs(x = NULL, y = "Mean performancce and 95% CI") +
        # remove x-axis coordiante extension
        scale_x_discrete(expand = c(0.05, 0)) +
        theme(panel.spacing = unit(0.1, "cm")) +
        theme(strip.text = element_text(size = 6)) +
        theme(axis.text = element_text(size = 6)) +
        # remove legend title
        theme(axis.title = element_text(size = 6)) +
        theme(legend.position = "right") +
        theme(legend.text = element_text(size = 6), legend.title = element_blank()) +
        theme(legend.key = element_rect(fill = "transparent", color = "transparent")) +
        theme(
                legend.box.background = element_rect(fill = "transparent", color = "transparent"),
                legend.box.spacing = unit(0, "pt"), # The spacing between the plotting area and the legend box (unit)
                legend.margin = margin(0, 0, 0, 0),
                plot.margin = margin(0, 0, 0, 0)
        )



ggsave(main_fig_between_cnns,
        p_main,
        width = 75,
        height = 60,
        unit = plot_unit,
        create.dir = TRUE
)

message(main_fig_between_cnns)

# save the plot to rds file
plot_file_rds <- gsub(".pdf", ".rds", main_fig_between_cnns)
saveRDS(p_main, plot_file_rds)

# hide legend
p_main_no_legend <- p_main + theme(legend.position = "none")
main_fig_between_cnns_no_legend <- gsub(".pdf", "_no_legend.pdf", main_fig_between_cnns)
ggsave(main_fig_between_cnns_no_legend,
        p_main_no_legend,
        width = 50,
        height = 60,
        unit = plot_unit,
        create.dir = TRUE
)
message(main_fig_between_cnns_no_legend)







################################################################################
# plot the boxplot comparison plot for each strata
################################################################################
ans <- xgb_cnn %>%
        filter(.metric %in% which_metrics_to_plot) %>%
        filter(ichorcna_strat == which_tf_strat) %>%
        filter(model_name %in% model_to_keep) %>%
        filter(feat %in% feat_chosen) %>%
        mutate(model_name_label = case_when(
                model_name == "xgboost" ~ "XGBoost",
                model_name == "resnet18" ~ "ResNet18",
                model_name == "lr" ~ "LR",
                TRUE ~ model_name
        ))

# add a col called '95% CI' to the data frame
ans <- ans %>%
        group_by(.metric, model_name_label) %>%
        mutate(
                .mean = Hmisc::smean.cl.boot(.estimate)["Mean"],
                .sd = Hmisc::smean.sd(.estimate)["SD"],
                .lower = Hmisc::smean.cl.boot(.estimate)["Lower"],
                .upper = Hmisc::smean.cl.boot(.estimate)["Upper"]
        )

# add a new col called 'anno' to the data frame

ans <- ans %>%
        group_by(.metric, model_name_label) %>%
        mutate(anno = paste("Mean: ", round(.mean, 3), "\n", "95% CI: ", round(.lower, 3), "-", round(.upper, 2), sep = ""))


ans$model_name_label <- factor(ans$model_name_label, levels = c("LR", "XGBoost", "ResNet18"))

# add col called 'model_name_label_anno' to the data frame
ans <- ans %>%
        group_by(.metric, model_name_label) %>%
        mutate(
                model_name_label_anno = paste(model_name_label, "\n", "Mean: ", round(.mean, 3), "\n", round(.lower, 3), "-", round(.upper, 2), sep = "")
        )

model_name_label_anno_vec <- ans$model_name_label_anno |>
        unique() |>
        as.character()

anno_data <- ans %>%
        group_by(.metric, model_name_label) %>%
        reframe(anno = paste("Mean: ", round(.mean, 3), "\n", round(.lower, 3), "-", round(.upper, 2), sep = "")) %>%
        distinct() %>%
        mutate(.estimate = 0.97) # 0.97 sets the position of the text annotation in each panel




# set the .estimate factor level to the same as which_metrics_to_plot
ans$.metric <- factor(ans$.metric, levels = which_metrics_to_plot)

p <- ggplot(ans, aes(x = model_name_label, y = .estimate, color = model_name_label)) +
        # stat_halfeye(
        #         # adjust bandwidth
        #         adjust = 0.5,
        #         # move to the right
        #         justification = -0.2,
        #         # remove the slub interval
        #         .width = 0,
        #         point_colour = NA
        # ) +
        # boxwhisker plot
        geom_violin(fill = "lightgrey", color = NA, width = 0.6) +
        geom_jitter(width = 0.2, height = 0, size = 0.5, alpha = 0.5) +
        geom_boxplot(outlier.shape = NA, color = "black", alpha = 0, width = 0.25, linewidth = 0.25) +
        # add annotation
        geom_text(
                data = anno_data,
                aes(label = anno),
                position = position_dodge(width = 0.9),
                size = 1.6,
                color = "black",
                vjust = 1.8
        ) +
        # add statistical annotation
        stat_compare_means(aes(x = model_name_label), method = "wilcox.test", label = "p.signif", ref.group = "LR", vjust = 0.8) +
        labs(
                y = y_axis_lab,
                # remove x lab
                x = NULL
        ) +
        theme_classic() +
        theme(
                axis.text.x = element_text(size = 5.5),
                axis.text.y = element_text(size = 5.5)
        ) +
        scale_color_manual(values = manual_colors) +
        facet_wrap(~ factor(.metric, levels = which_metrics_to_plot), scales = "free", nrow = 1, strip.position = "top") +
        # facet label color to lightgrey, size to 6
        theme(strip.background = element_rect(color = NA)) +
        # set strip panel color to lightgrey
        theme(strip.background = element_rect(fill = "lightgrey")) +
        # strip text to 5
        theme(strip.text = element_text(size = 5.5)) +
        # axis title to 5.5
        theme(axis.title = element_text(size = 5.5)) +
        # hide legend
        theme(legend.position = "none")

ggsave(plotfile, p, width = 18, height = 4.4, units = "cm", create.dir = TRUE, dpi = 300)

message("Saved plot to ", plotfile)



################################################################################
# plot AUROC, AUPRC, Sensitivity and Specificity
################################################################################
# filter the data
which_metrics_to_plot_in_optimal <- c("AUROC", "AUPRC", "Sensitivity", "Specificity")
ans <- xgb_cnn %>%
        filter(.metric %in% which_metrics_to_plot_in_optimal) %>%
        filter(ichorcna_strat == which_tf_strat) %>%
        filter(model_name %in% model_to_keep) %>%
        filter(feat %in% feat_chosen) %>%
        mutate(model_name_label = case_when(
                model_name == "xgboost" ~ "XGBoost",
                model_name == "resnet18" ~ "ResNet18",
                model_name == "lr" ~ "LR",
                TRUE ~ model_name
        ))

ans$.metric <- factor(ans$.metric, levels = which_metrics_to_plot_in_optimal)


# add a col called '95% CI' to the data frame
# ans <- ans %>%
#         group_by(.metric, model_name_label) %>%
#         mutate(
#                 .mean = mean(.estimate),
#                 .sd = sd(.estimate),
#                 .lower = .mean - 1.96 * sd(.estimate) / sqrt(n()),
#                 .upper = .mean + 1.96 * sd(.estimate) / sqrt(n())
#         )

ans <- ans %>%
        group_by(.metric, model_name_label) %>%
        mutate(
                .mean = Hmisc::smean.cl.boot(.estimate)["Mean"],
                .sd = Hmisc::smean.sd(.estimate)["SD"],
                .lower = Hmisc::smean.cl.boot(.estimate)["Lower"],
                .upper = Hmisc::smean.cl.boot(.estimate)["Upper"]
        )


# add a new col called '95% CI' to the data frame

ans <- ans %>%
        group_by(.metric, model_name_label) %>%
        mutate(anno = paste("Mean: ", round(.mean, 3), "\n", "95% CI: ", round(.lower, 3), "-", round(.upper, 2), sep = ""))


ans$model_name_label <- factor(ans$model_name_label, levels = c("LR", "XGBoost", "ResNet18"))

# add col called 'model_name_label_anno' to the data frame
ans <- ans %>%
        group_by(.metric, model_name_label) %>%
        mutate(
                model_name_label_anno = paste(model_name_label, "\n", "Mean: ", round(.mean, 3), "\n", round(.lower, 3), "-", round(.upper, 2), sep = "")
        )

model_name_label_anno_vec <- ans$model_name_label_anno |>
        unique() |>
        as.character()

anno_data <- ans %>%
        group_by(.metric, model_name_label) %>%
        reframe(anno = paste("Mean: ", round(.mean, 3), "\n", round(.lower, 3), "-", round(.upper, 2), sep = "")) %>%
        distinct() %>%
        mutate(.estimate = 0.97)




# set the .estimate factor level to the same as which_metrics_to_plot
ans$.metric <- factor(ans$.metric, levels = which_metrics_to_plot_in_optimal)

pvalue_vjust <- 0.8






if (which_tf_strat == "(0.1, 1]") {
        p <- ggplot(ans, aes(x = model_name_label, y = .estimate, color = model_name_label)) +
                # stat_halfeye(
                #         # adjust bandwidth
                #         adjust = 0.5,
                #         # move to the right
                #         justification = -0.2,
                #         # remove the slub interval
                #         .width = 0,
                #         point_colour = NA
                # ) +
                # boxwhisker plot
                geom_violin(fill = "lightgrey", color = NA, width = 0.9) +
                geom_jitter(width = 0.2, height = 0, size = 0.5, alpha = 0.5) +
                geom_boxplot(outlier.shape = NA, color = "black", alpha = 0, width = 0.35, linewidth = 0.25) +
                # add annotation
                geom_text(
                        data = anno_data,
                        aes(label = anno, group = model_name_label),
                        # position = position_dodge(width = 0.9),
                        size = 1.6,
                        vjust = 2.5,
                        color = "black"
                ) +
                # add statistical annotation
                stat_compare_means(aes(x = model_name_label),
                        method = "wilcox.test",
                        label = "p.signif",
                        ref.group = "LR",
                        # label.y = repeat(1.03, length(which_metrics_to_plot_in_optimal)),
                        label.y.npc = 0,
                        size = 3
                ) +
                labs(
                        y = y_axis_lab,
                        # remove x lab
                        x = NULL
                ) +
                theme_classic() +
                theme(
                        axis.text.x = element_text(size = 5.5),
                        axis.text.y = element_text(size = 5.5)
                ) +
                scale_color_manual(values = manual_colors) +
                facet_wrap(~ factor(.metric, levels = which_metrics_to_plot_in_optimal), scales = "fixed", nrow = 1, strip.position = "top") +
                # facet label color to lightgrey, size to 6
                theme(strip.background = element_rect(color = NA)) +
                # set strip panel color to lightgrey
                theme(strip.background = element_rect(fill = "lightgrey")) +
                # strip text to 5
                theme(strip.text = element_text(size = 5.5)) +
                # axis title to 5.5
                theme(axis.title = element_text(size = 5.5)) +
                # hide legend
                theme(legend.position = "none")
} else {
        p <- ggplot(ans, aes(x = model_name_label, y = .estimate, color = model_name_label)) +
                # stat_halfeye(
                #         # adjust bandwidth
                #         adjust = 0.5,
                #         # move to the right
                #         justification = -0.2,
                #         # remove the slub interval
                #         .width = 0,
                #         point_colour = NA
                # ) +
                # boxwhisker plot
                geom_violin(fill = "lightgrey", color = NA, width = 0.9) +
                geom_jitter(width = 0.2, height = 0, size = 0.5, alpha = 0.5) +
                geom_boxplot(outlier.shape = NA, color = "black", alpha = 0, width = 0.35, linewidth = 0.25) +
                # add annotation
                geom_text(
                        data = anno_data,
                        aes(label = anno),
                        position = position_dodge(width = 0.9),
                        size = 1.6,
                        color = "black",
                        vjust = 2.5
                ) +
                # add statistical annotation
                stat_compare_means(aes(x = model_name_label),
                        method = "wilcox.test",
                        label = "p.signif",
                        ref.group = "LR",
                        vjust = pvalue_vjust,
                        size = 3
                ) +
                labs(
                        y = y_axis_lab,
                        # remove x lab
                        x = NULL
                ) +
                theme_classic() +
                theme(
                        axis.text.x = element_text(size = 5.5),
                        axis.text.y = element_text(size = 5.5)
                ) +
                scale_color_manual(values = manual_colors) +
                facet_wrap(~ factor(.metric, levels = which_metrics_to_plot_in_optimal), scales = "free", nrow = 1, strip.position = "top") +
                # facet label color to lightgrey, size to 6
                theme(strip.background = element_rect(color = NA)) +
                # set strip panel color to lightgrey
                theme(strip.background = element_rect(fill = "lightgrey")) +
                # strip text to 5
                theme(strip.text = element_text(size = 5.5)) +
                # axis title to 5.5
                theme(axis.title = element_text(size = 5.5)) +
                # hide legend
                theme(legend.position = "none")
}



plotfile_new <- gsub(".pdf", "_sen_spe_optimal.pdf", plotfile)
ggsave(plotfile_new, p, width = 18, height = 4.4, units = "cm", create.dir = TRUE, dpi = 300)

message("Saved plot to ", plotfile_new)
