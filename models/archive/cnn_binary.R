
library(keras)
library(tensorflow)
library(tidyverse)
library(caret)
library(rsample)
library(reticulate)
library(MultiAssayExperiment)
library(EBImage)
options(warn = 1)

# parameters --------------------------------------------------------------
layers <- c("n_isize", "n_motif_s1_C", "n_motif_s1_T")
authors <- c("S.C. et al", "V.A. et al", "F.M. et al")
ichorcna_tf_cutoff_lower <- 0
ichorcna_tf_cutoff_upper <- 0.01
timepoints <- 1
mae_file <- "/mnt/scratchc/nrlab/wang04/ulyses/MAE3.rds"

outdir <- "/scratchc/nrlab/wang04/ulyses"

x_filename <- paste("x_clean", paste0(layers, collapse = "_"),"TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
y_filename <- paste("y_clean", paste0(layers, collapse = "_"),"TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
z_filename <- paste("z_clean", paste0(layers, collapse = "_"),"TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")

x_file <- file.path(outdir, x_filename)
y_file <- file.path(outdir, y_filename)
z_file <- file.path(outdir, z_filename)

plot_file_name <- file.path(outdir, "cnn_training_history_2.pdf")

# function ----------------------------------------------------------------

func_tibble_to_matrix <- function(x){
  # split the rowname column
  ans <- x %>%
    dplyr::mutate(bin = stringr::str_extract(rowname, "\\d+(p|q)_\\d+")) %>%
    dplyr::mutate(isize = stringr::str_extract(rowname, "isize_\\d+$")) %>%
    select(-rowname)
  # convert to matrix
  ans <- ans %>%
    pivot_wider(names_from = isize, values_from = value) %>%
    map_df(rev) %>%
    column_to_rownames("bin") %>% 
    as.matrix()
  return(ans)
}



func_split_by_assay <- function(x){
  ans <- x %>% 
    group_by(assay) 
  nms <- group_keys(ans)
  result <- ans %>%  
    group_split(.keep = FALSE) %>%
    purrr::set_names(nms[["assay"]])
  final <- purrr::map(result, func_tibble_to_matrix)
  return(final)
}


func_set_tensor_dimname <- function(x, dimname_to_set){
  dimnames(x) <- dimname_to_set
  return(x)
}


func_get_label <- function(bam_id, col_Data){
  ans <- col_Data[rownames(col_Data) == bam_id, ]$"cohort"
  return(ans)
}


func_get_patient_id <- function(bam_id, col_Data){
  ans <- col_Data[rownames(col_Data) == bam_id, ]$"patient_id"
  return(ans)
}

# select data ----------------------------------------------------------

# if warn = 2, it will treat warnings as errors

mae <- readRDS(mae_file)
author_boolean <- mae$author %in% authors
ichorcna_tf_boolean <- mae$ichorcna_tf < ichorcna_tf_cutoff_upper & mae$ichorcna_tf >= ichorcna_tf_cutoff_lower & !is.na(mae$ichorcna_tf)
timepoint_boolean <- mae$timepoint %in% timepoints | is.na(mae$timepoint)
healthy_boolean <- mae$cohort == "Healthy"

col_boolean <- (author_boolean & ichorcna_tf_boolean) | (author_boolean & healthy_boolean)
ans <- mae[, col_boolean, layers]
#ans_long <- longFormat(ans, colDataCols = c("cohort", "patient_id"))
ans_long <- longFormat(ans, colDataCols = "cohort")

# convert to matrix ------------------------------------------------------

# group by primary
ans_tibble<- ans_long %>%
  as_tibble() %>%
  select(-colname) %>%
  select(-cohort) %>%
  group_by(primary)

ans_primary <- ans_tibble %>%
  group_split(.keep = FALSE) %>%
  set_names(group_keys(ans_tibble)[["primary"]])


# use future map to speed up
ans_primary_assay <- purrr::map(ans_primary, func_split_by_assay)
ans_primary_assay_tensor <- purrr::map(ans_primary_assay, EBImage::combine)

dimnames_to_set <- list(ans_primary_assay_tensor[[1]] %>% 
                            dimnames() %>% pluck(1), 
                        ans_primary_assay_tensor[[1]] %>% 
                          dimnames() %>% 
                          pluck(2),
                        names(ans_primary_assay[[1]])
                        )

tensors <-  purrr::map(ans_primary_assay_tensor,
                      .f = func_set_tensor_dimname,
                      dimname_to_set = dimnames_to_set)

# get labels
bam_ids <- names(tensors)
mae_coldata <- colData(mae)
labels <- purrr::map(bam_ids, .f = func_get_label, col_Data  = mae_coldata)

pids <- purrr::map(bam_ids, .f = func_get_patient_id, col_Data  = mae_coldata)

# convert to tensor ------------------------------------------------------
# shuffle the data randomly
shuffle_index <- sample(seq_along(labels))
tensors <- tensors[shuffle_index]
labels <- labels[shuffle_index]
pids <- pids[shuffle_index]

x <- abind::abind(tensors, along = 0)

y <- unlist(labels) %>%
  as_tibble() %>%
  mutate(cohort = case_when(
    value == "Healthy" ~ 0L,
    TRUE ~ 1L
  )) %>%
  purrr::pluck("cohort")

z <- unlist(pids)

# save x and y as npy files using numpy
np <- reticulate::import("numpy")

np$save(x_file, x)
np$save(y_file, y)
np$save(z_file, z)

# copy file to /Users/wang04, with replace
file.copy(x_file, "/Users/wang04", overwrite = TRUE)
file.copy(y_file, "/Users/wang04", overwrite = TRUE)
file.copy(z_file, "/Users/wang04", overwrite = TRUE)

# run the model ----------------------------------------------------------
# run model cnn_binary_cross_validation.py using python

# reticulate run python script via singularity container
reticulate::use_condaenv("ulyses", required = TRUE)
singularity run --nv /mnt/scratchc/nrlab/wang04/ulyses/ulyses.sif  /mnt/scratchc/nrlab/wang04/ulyses/cnn_binary_cross_validation.py --x_file /mnt/scratchc/nrlab/wang04/ulyses/x_clean.npy --y_file /mnt/scratchc/nrlab/wang04/ulyses/y_clean.npy --plot_file_name /mnt/scratchc/nrlab/wang04/ulyses/cnn_training_history_2.pdf --epochs 100 --batch_size 32 --n_splits 5 --n_repeats 10 --n_jobs 10 
reticulate::py_run_file("cnn_binary_cross_validation.py")

# get results and plot in R ----------------------------------------------
## access the python variable

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














# split samples -----------------------------------------------------------

# train_ids <- caret::createDataPartition(y, p = 0.67)
# split_result <- rsample::initial_split(tibble("cohort" = y), strata="cohort", prop = 0.67)
# train_ids <- split_result$in_id

sklearn <- reticulate::import("sklearn.model_selection")

split_result <- sklearn$train_test_split(x, y, test_size=0.2, shuffle=TRUE, stratify = y) %>%
  purrr::set_names(c("x_train_validation", "x_test", "y_train_validation", "y_test"))


x_train_validation <- split_result$x_train_validation
y_train_validation <- split_result$y_train_validation


split_train_validation_result <- sklearn$train_test_split(x_train_validation, y_train_validation, test_size=0.45, shuffle=TRUE, stratify = y_train_validation) %>%
  purrr::set_names(c("x_train", "x_validation", "y_train", "y_validation"))

x_train <- split_train_validation_result$x_train
x_validation <- split_train_validation_result$x_validation

y_train <- split_train_validation_result$y_train
y_validation <- split_train_validation_result$y_validation

x_test <- split_result$x_test
y_test <- split_result$y_test



# set class weight --------------------------------------------------------
class_weight <-class_weights(y_train) 
input_shape <- dim(tensors[[1]])

# build model -------------------------------------------------------------

# convolutional layer
model <- keras_model_sequential() %>% 
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu", 
                input_shape = input_shape) %>% 
  layer_max_pooling_2d(pool_size = c(3,2)) %>% 
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu") %>% 
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu") %>% 
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu") %>% 
  layer_max_pooling_2d(pool_size = c(2,2)) 

model %>% 
  layer_flatten() %>% 
  layer_dense(units = 32, activation = "relu") %>% 
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  optimizer = optimizer_adam(learning_rate = 0.0008),
  loss = "binary_crossentropy",
  metrics = list("accuracy", 
                 "AUC", 
                 metric_sensitivity_at_specificity(specificity = 0.95, name = "sen_at_0.95spec"), 
                 metric_sensitivity_at_specificity(specificity = 0.98, name = "sen_at_0.98spec"), 
                 metric_precision(), 
                 metric_recall())
)


# train the model
history_cnn <- model %>% 
  fit(
    x = x_train, 
    y = y_train,
    epochs = 40,
    batch_size = 60,
    verbose = 1,
    callbacks = callback_tensorboard(),
    validation_data = list(x_validation, y_validation),
    shuffle = TRUE,
    class_weight = class_weight
  )


# evaluate the model

model %>% evaluate(x_test, y_test)

pred <- model %>% predict(x_test)

#plot the history
# pdf(file = plot_file_name, width=8, height=12)

p <- plot(history_cnn) 
ggsave(filename = plot_file_name, plot = p,  width = 7 , height = 10)


