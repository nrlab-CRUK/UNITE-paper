library(tidyverse)
#https://stackoverflow.com/questions/52467915/plotting-mean-roc-curve-for-multiple-roc-curves-r
input_path <- "/Users/wang04/ulyses/models/bill/roc_curve_efficient_net_s_tf_0_0.01/"
roc_file_pattern <- "_y_test_y_pred_z_test"

fl <- list.files(input_path, pattern = roc_file_pattern, full.names = TRUE, recursive = TRUE)


ans <- lapply(fl, read_csv )
names(ans) <- as.vector(fl)
ans2 <- bind_rows(ans, .id = "file")


#-----

library(cutpointr)
library(tidyverse)
mean_roc <- function(data, cutoffs = seq(from = 0, to = 1, by = 0.01)) {
  map_df(cutoffs, function(cp) {
    out <- cutpointr(data = ans2, x = y_pred, class = y_test,
                     subgroup = file, method = oc_manual, cutpoint = cp,
                     pos_class = 1, direction = ">=")
    data.frame(cutoff = cp, 
               sensitivity = mean(out$sensitivity),
               specificity = mean(out$specificity))
  })
}

mr <- mean_roc(ans2)





#####metrics
library(cutpointr)
library(tidyverse)


metrics_file_pattern <- "_metrics.csv"
metrics_fl <- list.files(input_path, pattern = metrics_file_pattern, full.names = TRUE, recursive = TRUE)

suppressWarnings(
  metric_fl_read <- lapply(metrics_fl, read_csv) %>%
    bind_rows() %>%
    # remove duplicated rows
    distinct()
  
)

sen_99_spec <- mean(metric_fl_read$test_sen_99spe) |> 
  round(2)

auroc_mean <- mean(metric_fl_read$test_auroc) |> 
  round(2)
auroc_sd <- sd(metric_fl_read$test_auroc) |>
  round(2)

# create text for mean and sd, mean +- sd
auroc_text_ <- paste0("EfficientNet AUROC: ", auroc_mean, "±", auroc_sd)



# plot roc curve

p <- cutpointr(data = ans2, 
               x = y_pred, class = y_test, subgroup = file,
               pos_class = 1, direction = ">=") %>% 
  plot_roc(display_cutpoint = FALSE, type = "line") +
  # add plot title 
  ggtitle("ROC Curve") +
  theme_classic() +
  # remove legend position
  theme(legend.position = "none") +
  geom_line(data = mr, 
            mapping = aes(x = 1 - specificity, y = sensitivity), 
            color = "black", 
            linewidth = 1 ) +
  # add dashed line to as diagonal
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") 


p_mean <-   ggplot() +
  geom_line(data = mr, 
            mapping = aes(x = 1 - specificity, y = sensitivity), 
            color = "black", 
            linewidth = 1 ) +
  geom_vline(xintercept = 0.01, linetype = "dashed") +
  theme_classic() +
  coord_cartesian(xlim = c(0, 0.05))


# add text to plot p
p_anno <- p + annotate("text", x = 0.75, y = 0.25, label = auroc_text_, size = 3)

p_mean_anno <- p_mean


# save plot height. = 5cm width = 5cm, save as pdf
ggsave("/Users/wang04/ulyses/models/bill/roc_curve_efficient_net_s_tf_0_0.01/roc_curve_plot_all.pdf", 
       plot = p_anno, 
       height = 9, width = 9, units = "cm")

ggsave("/Users/wang04/ulyses/models/bill/roc_curve_efficient_net_s_tf_0_0.01/roc_curve_plot_zoom.pdf", 
       plot = p_mean_anno, 
       height = 5, width = 5, units = "cm")




