suppressMessages(library(keras))
suppressMessages(library(tensorflow))
suppressMessages(library(tidyverse))
suppressMessages(library(caret))
suppressMessages(library(rsample))
suppressMessages(library(reticulate))
suppressMessages(library(MultiAssayExperiment))
suppressMessages(library(EBImage))
suppressMessages(library(bettermc))
suppressMessages(library(furrr))
options(warn = 1)


# function ----------------------------------------------------------------

func_tibble_to_matrix <- function(x) {
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



func_split_by_assay <- function(x) {
  ans <- x %>%
    group_by(assay)
  nms <- group_keys(ans)
  result <- ans %>%
    group_split(.keep = FALSE) %>%
    purrr::set_names(nms[["assay"]])
  final <- purrr::map(result, func_tibble_to_matrix)
  return(final)
}


func_set_tensor_dimname <- function(x, dimname_to_set) {
  # if x only has two dimensions, then add a 3rd dimension
  if (length(dim(x)) == 2) {
    message("Adding a 3rd dimension to the tensor.")
    x <- array(x, dim = c(dim(x), 1))
  }

  dimnames(x) <- dimname_to_set
  return(x)
}


func_get_label <- function(bam_id, col_Data) {
  ans <- col_Data[rownames(col_Data) == bam_id, ]$"cohort"
  return(ans)
}


func_get_patient_id <- function(bam_id, col_Data) {
  ans <- col_Data[rownames(col_Data) == bam_id, ]$"patient_id"
  return(ans)
}



tensor_maker <- function(layers = c("n_isize", "n_motif_smono1_C", "n_motif_smono1_T"),
                         authors = c("V.A. et al", "F.M. et al"),
                         ichorcna_tf_cutoff_lower,
                         ichorcna_tf_cutoff_upper,
                         timepoints = NULL,
                         cohorts = NULL,
                         sample_type = c("plasma"),
                         mae_file,
                         seed = 0,
                         final_test_hyper = NULL,
                         outdir) {
  # libs ------------------------------------------------------------------
  suppressMessages(library(keras))
  suppressMessages(library(tensorflow))
  suppressMessages(library(tidyverse))
  suppressMessages(library(caret))
  suppressMessages(library(rsample))
  suppressMessages(library(reticulate))
  suppressMessages(library(MultiAssayExperiment))
  suppressMessages(library(EBImage))
  suppressMessages(library(bettermc))
  options(warn = 1)


  # set train filenames --------------------------------------------------------
  x_filename <- paste("x_clean", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  y_filename <- paste("y_clean", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  z_filename <- paste("z_clean", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  bam_id_filename <- paste("bam_id_clean", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  x_file <- file.path(outdir, x_filename)
  y_file <- file.path(outdir, y_filename)
  z_file <- file.path(outdir, z_filename)
  bam_id_file <- file.path(outdir, bam_id_filename)

  # set test filenames --------------------------------------------------------
  x_test_filename <- paste("x_clean_independent_test", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  y_test_filename <- paste("y_clean_independent_test", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  z_test_filename <- paste("z_clean_independent_test", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  bam_id_test_filename <- paste("bam_id_clean_independent_test", "layers", "TF", ichorcna_tf_cutoff_lower, ichorcna_tf_cutoff_upper, ".npy", sep = "_")
  x_test_file <- file.path(outdir, x_test_filename)
  y_test_file <- file.path(outdir, y_test_filename)
  z_test_file <- file.path(outdir, z_test_filename)
  bam_id_test_file <- file.path(outdir, bam_id_test_filename)


  plot_file_name <- file.path(outdir, "cnn_training_history_2.pdf")

  # read data ------------------------------------------------------------
  message("reading MAE data.")
  mae <- readRDS(mae_file)
  # filter data ----------------------------------------------------------
  sample_type_boolean <- mae$sample_type %in% sample_type
  author_boolean <- mae$author %in% authors
  # if ichorcna_tf_cutoff_upper=0 and ichorcna_tf_cutoff_upper=1, then use all values
  if (ichorcna_tf_cutoff_upper == 1 & ichorcna_tf_cutoff_lower == 0) {
    ichorcna_tf_boolean <-
      (mae$ichorcna_tf <= ichorcna_tf_cutoff_upper) &
        (mae$ichorcna_tf >= ichorcna_tf_cutoff_lower) &
        (!is.na(mae$ichorcna_tf))
    # if ichorcna_tf_cutoff_upper=0.01 and ichorcna_tf_cutoff_upper=0, then
  } else if (ichorcna_tf_cutoff_upper == 0.03 & ichorcna_tf_cutoff_lower == 0) {
    ichorcna_tf_boolean <-
      (mae$ichorcna_tf <= ichorcna_tf_cutoff_upper) &
        (mae$ichorcna_tf >= ichorcna_tf_cutoff_lower) &
        (!is.na(mae$ichorcna_tf))
  } else {
    ichorcna_tf_boolean <-
      (mae$ichorcna_tf <= ichorcna_tf_cutoff_upper) &
        (mae$ichorcna_tf > ichorcna_tf_cutoff_lower) &
        (!is.na(mae$ichorcna_tf))
  }

  # if timepoints is NULL, then use all timepoints
  if (is.null(timepoints)) {
    timepoint_boolean <- TRUE
  } else {
    timepoint_boolean <- mae$timepoint %in% timepoints
  }

  # if cohort is NULL, then use all cohorts
  if (is.null(cohorts)) {
    cohort_boolean <- TRUE
  } else {
    cohort_boolean <- mae$cohort %in% cohorts
  }

  healthy_boolean <- mae$cohort == "Healthy"

  outlier_boolean <- !rownames(colData(mae)) %in% outliers_hyper

  col_boolean <- (author_boolean & outlier_boolean & ichorcna_tf_boolean & timepoint_boolean & sample_type_boolean & cohort_boolean) |
    (author_boolean & outlier_boolean & healthy_boolean & timepoint_boolean & sample_type_boolean)

  message("filtering MAE data.")
  ans <- mae[, col_boolean, layers]
  message("Converting to long format.")
  ans_long <- longFormat(ans, colDataCols = "cohort")

  # convert to matrix ------------------------------------------------------

  # group by primary
  ans_tibble <- ans_long %>%
    as_tibble() %>%
    select(-colname) %>%
    select(-cohort) %>%
    group_by(primary)

  ans_primary <- ans_tibble %>%
    group_split(.keep = FALSE) %>%
    set_names(group_keys(ans_tibble)[["primary"]])


  # TODO use future map to speed up
  ans_primary_assay <- purrr::map(ans_primary, func_split_by_assay)
  ans_primary_assay_tensor <- purrr::map(ans_primary_assay, EBImage::combine)

  dimnames_to_set <- list(
    ans_primary_assay_tensor[[1]] %>%
      dimnames() %>% pluck(1),
    ans_primary_assay_tensor[[1]] %>%
      dimnames() %>%
      pluck(2),
    names(ans_primary_assay[[1]])
  )


  # if the 3rd dimension only contains one element, then remove it
  # if (length(dimnames_to_set[[3]]) == 1){
  #  dimnames_to_set <- dimnames_to_set[1:2]
  # }

  message("Setting tensor dimnames.")

  tensors <- purrr::map(ans_primary_assay_tensor,
    .f = func_set_tensor_dimname,
    dimname_to_set = dimnames_to_set
  )

  # get labels
  bam_ids <- names(tensors)
  mae_coldata <- colData(mae)
  labels <- purrr::map(bam_ids, .f = func_get_label, col_Data = mae_coldata)

  pids <- purrr::map(bam_ids, .f = func_get_patient_id, col_Data = mae_coldata)

  # convert to tensor ------------------------------------------------------
  # shuffle the data randomly
  # set seed
  set.seed(seed)
  shuffle_index <- sample(seq_along(labels))
  tensors <- tensors[shuffle_index]
  labels <- labels[shuffle_index]
  pids <- pids[shuffle_index]
  bam_ids <- bam_ids[shuffle_index]

  x <- abind::abind(tensors, along = 0)

  # if x only has 3 dimensions, then add a 4th dimension
  if (length(dim(x)) == 3) {
    x <- array(x, dim = c(dim(x), 1))
  }

  # issue a warning if their are more than 4 dimensions
  if (length(dim(x)) > 4) {
    warning("The tensor has more than 4 dimensions. Check the reason. usually 1d is n samples ")
  }

  y <- unlist(labels) %>%
    as_tibble() %>%
    mutate(cohort = case_when(
      value == "Healthy" ~ 0L,
      TRUE ~ 1L
    )) %>%
    purrr::pluck("cohort")

  z <- unlist(pids)

  # split into train and test

  load(final_test_hyper)

  # train set
  final_train_index <- which(!bam_ids %in% final_test_primary_vec)
  x_train <- x[final_train_index, , , ]
  y_train <- y[final_train_index]
  z_train <- z[final_train_index]
  bam_id_train <- bam_ids[final_train_index]

  # test set
  final_test_index <- which(bam_ids %in% final_test_primary_vec)
  x_test <- x[final_test_index, , , ]
  y_test <- y[final_test_index]
  z_test <- z[final_test_index]
  bam_id_test <- bam_ids[final_test_index]



  # save train files ------------------------------------------------------------

  use_condaenv("R4_2")
  np <- reticulate::import("numpy")
  # convert bam_id_train to a tibble
  print("saving train files")
  print("----------------")

  np$save(x_file, x_train)
  message("saved x to ", x_file)
  np$save(y_file, y_train)
  message("saved y to ", y_file)
  np$save(z_file, z_train)
  message("saved z to ", z_file)
  np$save(bam_id_file, bam_id_train)
  message("saved bam_id to ", bam_id_file)

  message("Done!")

  # save test files ------------------------------------------------------------
  print("saving test files")
  print("----------------")

  np$save(x_test_file, x_test)
  message("saved x to ", x_test_file)
  np$save(y_test_file, y_test)
  message("saved y to ", y_test_file)
  np$save(z_test_file, z_test)
  message("saved z to ", z_test_file)
  np$save(bam_id_test_file, bam_id_test)
  message("saved bam_id to ", bam_id_test_file)

  message("Done!")
}
