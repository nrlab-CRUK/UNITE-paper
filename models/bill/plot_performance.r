library(tidyverse)
library(patchwork)
input_path <- "/home/nrlab/wang04/ulyses/models/bill"
fl <- list.files(path = input_path, pattern = "performance_efficient_net_s_tf.*\\.csv", full.names = TRUE)

fl <- setNames(fl, c("0-0.01", "all", "0.01-0.03", "0.03-0.1", "0.1-0.2", "0.2-1"))

ans <- lapply(fl, readr::read_csv)

tb <- bind_rows(ans, .id = "tf")

# plot test_auroc as box-whisker plot, x is tf and y is test_auroc
p_auroc <- tb %>%
        ggplot(aes(x = tf, y = test_auroc)) +
        geom_boxplot(outlier.shape = NA) +
        # show median value beside the box
        stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)), hjust = -1) +
        geom_jitter(width = 0.2, height = 0.1, alpha = 0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "EfficientNet-S", x = "Transfer Factor", y = "Test AUROC") +
        scale_y_continuous(limits = c(0.7, 1))
p_acc <- tb %>%
        ggplot(aes(x = tf, y = test_acc)) +
        geom_boxplot(outlier.shape = NA) +
        # show median value beside the box
        stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)), hjust = -1) +
        geom_jitter(width = 0.2, height = 0.1, alpha = 0.5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "EfficientNet-S", x = "Transfer Factor", y = "Test Accuracy") +
        scale_y_continuous(limits = c(0.7, 1))

p <- p_auroc / p_acc

p
