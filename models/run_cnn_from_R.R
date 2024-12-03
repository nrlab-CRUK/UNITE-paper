
#############################################################################
# run model cnn_binary_cross_validation.py using python
#############################################################################

# reticulate run python script via singularity container
reticulate::use_condaenv("ulyses", required = TRUE)
singularity run --nv /mnt/scratchc/nrlab/wang04/ulyses/ulyses.sif  /mnt/scratchc/nrlab/wang04/ulyses/cnn_binary_cross_validation.py --x_file /mnt/scratchc/nrlab/wang04/ulyses/x_clean.npy --y_file /mnt/scratchc/nrlab/wang04/ulyses/y_clean.npy --plot_file_name /mnt/scratchc/nrlab/wang04/ulyses/cnn_training_history_2.pdf --epochs 100 --batch_size 32 --n_splits 5 --n_repeats 10 --n_jobs 10 
reticulate::py_run_file("cnn_binary_cross_validation.py")

# get results and plot in R ----------------------------------------------
cv_df <- py$cv_df

cv_df <- read_csv(file = "/home/wang04/cv_results_5fold_10repeats_n_isize.csv" )
cv_df <- read_csv(file = "/home/wang04/cv_results_5fold_10repeats.csv" )

# visualize the results


cv_long <- cv_df %>%
  pivot_longer(colnames(cv_df))


p <- ggplot(cv_long, aes(name, value, fill = name)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(size = 1,width = 0.4, show.legend = FALSE) +
  stat_summary(aes(group = name, label = round(after_stat(y),3)), fun = mean, geom = "text", color = "black", hjust = -2.4) +
  labs(x = "Score", y = "Value") +
  theme_classic() +
  theme(text = element_text(size = 20)) 

viz_filename <- "cnn_binary_cross_validation_5cv_10repeats_n_isize.pdf"

viz_filename <- "cnn_binary_cross_validation_5cv_10repeats.pdf"
ggsave(filename = viz_filename, plot = p, width = 15, height = 8, dpi = 300)

message("ploted: ", viz_filename)

