library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)


# parameters
tmad_vs_ichor_file <- "/home/nrlab/wang04/ulyses/plot_mae/tmad_vs_ichorcna_tf.pdf"
triple_plot_file <- "/home/nrlab/wang04/ulyses/plot_mae/predicted_VS_exprected.pdf"

source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
mae <- readRDS(mae_file_latest)

col_data <- colData(mae)
tmad_ichor <- tibble(tmad = col_data$tmad, ichorcna_tf = col_data$ichorcna_tf)



################################################################################
# plot cor between tMAD and ichorcna_tf as a scatter plot with trend line and R2
################################################################################

x_breaks <- c(0, 0.03, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
y_breaks <- c(0, 0.015, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

tmad_ichor <- tmad_ichor %>%
  dplyr::filter_all(all_vars(!is.na(.)))


p_tmad_ichor <- ggplot(tmad_ichor, 
              aes(y= tmad, 
                  x = ichorcna_tf)) +
  geom_point(color = "black", 
              fill = "#e2dede", 
              shape = 21, 
              size = 1 ) +
  # add a trend line
  geom_smooth(method = "lm", 
              se = TRUE, 
              color = "black", 
              linetype = "dashed") +
  # add R2 to the plot
  stat_cor(method = "pearson", 
          label.x.npc = 0.1,
          label.y = 0.9,
           label.sep = " , ", 
           size = 3.5) +
  # add text annoation about the number of samples
  annotate("text", 
           x = 0.2, 
           y = 0.85, 
           label = paste("( n = ", nrow(all_tmad_ichor), ")", sep = "")) +

  # add a line at y = 0.015
  geom_hline(yintercept = 0.015, 
             linetype = "dashed", 
             color = "red") +
  # add a line at x = 0.03
  geom_vline(xintercept = 0.03, 
             linetype = "dashed", 
             color = "red") +
  
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = x_breaks,
                     labels = as.character(x_breaks),
                     limits = c(0.0, 1.0)) +
  scale_y_continuous(breaks = y_breaks,
                     labels = as.character(y_breaks),
                     limits = c(0.0, 1.0)) +
  labs(y = "tMAD", x = "iChorCNA TF") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # make axis text smaller
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  theme(strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"))

ggsave(tmad_vs_ichor_file, p_tmad_ichor, width = 5, height = 5.1)

message(tmad_vs_ichor_file)


################################################################################
# plot triple_filtered$tmad, triple_filtered$clinical_tf, triple_filtered$ichorcna_tf in a scatter plot
################################################################################


triple_filtered <- tibble(tMAD = col_data$tmad, 
			'ichorCNA TF' = col_data$ichorcna_tf, 
			'Expected TF' = col_data$clinical_tf) %>%
  dplyr::filter_all(all_vars(!is.na(.)))


triple_filtered_long <- triple_filtered %>%
  pivot_longer(cols = c("tMAD", "ichorCNA TF"), 
               names_to = "metrics", 
               values_to = "predicted") %>%
  mutate(variable = factor(metrics, levels = c('tMAD', 'ichorCNA TF')))


x_breaks <- c(0, 0.015, 0.03, 0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
y_breaks <- c(0,  0.1,  0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)


p <- ggplot(triple_filtered_long, 
              aes(x = predicted, 
                  y = `Expected TF`, 
                  color = metrics)) +
  geom_point(color = "black", 
              fill = "#e2dede", 
              shape = 21, 
              size = 2 ) +
  # add a trend line
  geom_smooth(method = "lm", 
              se = TRUE, 
              color = "black", 
              linetype = "dashed") +
  # add R2 to the plot
  stat_cor(method = "pearson", 
          label.x.npc = 0.2,
          label.y = 0.9,
           label.sep = " , ", 
           size = 3.5) +
  # add text annoation about the number of samples
  annotate("text", 
           x = 0.12, 
           y = 0.85, 
           label = paste("( n = ", nrow(triple_filtered), ")", sep = "")) +
  # add vertical line at x = 0.03
  geom_vline(xintercept = 0.03, 
	     linetype = "dashed", 
	     color = "red") +
  # add vertical line at x = 0.015
  geom_vline(xintercept = 0.015, 
	     linetype = "dashed", 
	     color = "red") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = x_breaks,
                     labels = as.character(x_breaks),
                     limits = c(0.0, 0.41)) +
  scale_y_continuous(breaks = y_breaks,
                     labels = as.character(y_breaks),
                     limits = c(0.0, 1.0)) +
  facet_wrap(~variable) +
  labs(x = "Score", y = "Expected TF") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_rect(fill = "lightgrey", colour = "lightgrey"))



ggsave(triple_plot_file, p, width = 11, height = 4)
message(triple_plot_file)

################################################################################
# alluvial plot
################################################################################
