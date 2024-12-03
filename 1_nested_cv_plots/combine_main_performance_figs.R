library(patchwork)
library(tidyverse)

xgb_benchmark_file <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/xgboost_ranking_plot_main_fig.rds"
cnn_benchmark_file <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/resnet18/C4/cnn_ranking_plot_main_fig.rds"
xgb_vs_cnn_file <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_and_cnn/xgb_vs_cnn_main_fig.rds"
cnn_internal_file <- "/home/nrlab/wang04/ulyses/1_nested_cv_plots/direct_comp_xgb_and_cnn/cnn_models_internal_comparison.rds"

xgb_benchmark <- readRDS(xgb_benchmark_file)
xgb_benchark <- xgb_benchmark 
cnn_benchmark <- readRDS(cnn_benchmark_file)
cnn_benchmark <- cnn_benchmark 
xgb_vs_cnn <- readRDS(xgb_vs_cnn_file) 
cnn_internal <- readRDS(cnn_internal_file) 

p <- xgb_benchmark +
        cnn_benchmark +
        xgb_vs_cnn +
        cnn_internal +
        plot_layout(ncol = 2, guides = "keep")
p2 <- p &
        # set the plot margins to default
        theme(plot.margin = margin(1, 1, 1, 1, "mm")) &
        # set axis title size to 6
        theme(axis.title = element_text(size = 6.5))

# save the plot
ggsave("/home/nrlab/wang04/ulyses/1_nested_cv_plots/main_performance_figs.pdf", p2,
        width = 18, height = 10, units = "cm", dpi = 300
)
message("Saved the main performance figures to /home/nrlab/wang04/ulyses/1_nested_cv_plots/main_performance_figs.pdf
")
