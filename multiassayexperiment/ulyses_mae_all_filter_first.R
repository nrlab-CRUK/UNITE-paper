# libs --------------------------------------------------------------------
library(tidyverse)
library(MultiAssayExperiment)
library(rslurm)

# parameters --------------------------------------------------------------

data_path <- "/mnt/scratchc/nrlab/wang04"
MAE_output_path <- "/scratchc/nrlab/wang04/build_mae"
re_read_ulyses_obj <- TRUE
save_ulyses_obj <- FALSE
meta_file <- "/home/nrlab/wang04/ulyses/meta_data/colData_batch1_2.rds"
meta_file_csv <- "/home/nrlab/wang04/ulyses/meta_data/colData_batch1_2.csv"
if (TRUE) {
  message("Reading in ", meta_file_csv, " ...")
  meta_data_tmp <- read.csv(meta_file_csv)
  saveRDS(meta_data_tmp, meta_file)
  message("Meta data saved to ", meta_file, " ...")
}

authors <- c("all")

# make author label for file name
author_label <- authors |>
  paste(collapse = "_") |>
  stringr::str_replace_all("[[:punct:]]", "") |>
  stringr::str_replace_all(" ", "_") |>
  # remove the _et_al
  stringr::str_replace_all("_et_al", "")

ulyses_obj_rds_stored <- file.path("/mnt/scratchc/nrlab/wang04/build_mae", paste0(basename(meta_file), "_ulyses_obj_", author_label, ".rds"))
mae_filename <- file.path(MAE_output_path, paste0(basename(meta_file), "_author_", author_label, "_MAE.rds"))

# read in the meta file
ulyses_meta <- readRDS(meta_file) %>%
  dplyr::filter(!is.na(bam_id)) %>%
  dplyr::filter(bam_id != "SLX-10991.NA") %>%
  dplyr::filter(sample_type == "plasma") %>%
  dplyr::filter(timepoint == 1) %>%
  # remove when ichorcna_tf is NA
  dplyr::filter(!is.na(ichorcna_tf)) %>%
  # remove when author is "Stephen Cristiano"
  dplyr::filter(author != "Stephen Cristiano")
# TODO: add filter for seq_depth

# only keep author %in% authors
# skip this step if authors is "all"
if (!identical(authors, "all")) {
  ulyses_meta <- ulyses_meta %>%
    dplyr::filter(author %in% authors)
}

# functions ---------------------------------------------------------------

readRDS_hw <- function(x) {
  message("Reading in ", x, " ...")
  obj <- readRDS(x)
}

get_bam_id <- function(x) {
  if (stringr::str_detect(x, "EE\\d+\\.hg")) {
    bam_id <- stringr::str_extract(basename(x), "(^EE\\d+)\\.hg\\d\\d", group = 1)
  } else {
    bam_id <- stringr::str_extract(basename(x), "(^(SLX|DL)\\-\\d+\\.(\\w+\\d+(\\w+)?(\\-\\w+\\d+(\\w+)?)?)|(\\w+))\\.", group = 1)
  }
  return(bam_id)
}

# readin rds files --------------------------------------------------------

if (re_read_ulyses_obj) {
  file_list <- list.files(path = data_path, pattern = "ulyses_obj.rds", full.names = TRUE, recursive = TRUE) %>%
    as.list()

  # remove files containing "archive" or "Archive" or "test" or "debug" in the path
  file_list <- file_list[!grepl("archived|archive|Archive|Archived|test|Test|debug|Debug", file_list)]
  # remove files containing "csf" or "urine" in the path
  file_list <- file_list[!grepl("csf|urine", file_list)]

  # filter out the files that are not in the meta file bam_id col
  bam_id_all <- lapply(file_list, get_bam_id)
  file_list_filter <- file_list[which(bam_id_all %in% ulyses_meta$bam_id)]
  message("Reading in ", length(file_list_filter), " files ...")
  file_list_rds <- bettermc::mclapply(file_list_filter, readRDS, mc.cores = 16)
} else {
  file_list_rds <- readRDS(file = ulyses_obj_rds_stored)
}

if (save_ulyses_obj) {
  message("Saving ulyses obj ...")
  saveRDS(object = file_list_rds, file = ulyses_obj_rds_stored)
  message("Ulyses obj saved to ", ulyses_obj_rds_stored, " ...")
}

# functions ---------------------------------------------------------------

get_bam_id <- function(x) {
  if (stringr::str_detect(x, "EE\\d+\\.hg")) {
    bam_id <- stringr::str_extract(basename(x), "(^EE\\d+)\\.hg\\d\\d", group = 1)
  } else {
    bam_id <- stringr::str_extract(basename(x), "(^(SLX|DL)\\-\\d+\\.(\\w+\\d+(\\w+)?(\\-\\w+\\d+(\\w+)?)?)|(\\w+))\\.", group = 1)
  }
  return(bam_id)
}

lookup_primary <- function(x, input_col, output_col, meta_file = ulyses_meta) {
  meta_file[which(meta_file[input_col] == x), ][[output_col]]
}

lookup_primary_boolean <- function(x, input_col, meta_file = ulyses_meta) {
  if (x %in% meta_file[[input_col]]) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


tensor_helper <- function(x) {
  dplyr::mutate(x, bin_isize = paste(Var1, Var2, sep = "_")) |>
    dplyr::select(bin_isize, value)
}


tensor_helper_merge <- function(x, tensor_assay, tensor_coldata) {
  tensor_helper_get_layer <- function(input, layer) {
    if (layer %in% names(input)) {
      message("Handling", "   ", layer, "...")
      ans <- input[[layer]] |>
        column_to_rownames("bin_isize") |>
        as.matrix()
      return(ans)
    } else {
      ans <- input[["n_isize"]] |>
        column_to_rownames("bin_isize") |>
        as.matrix()
      ans[, 1] <- NA
      return(ans)
    }
  }


  final <- map(tensor_assay, .f = tensor_helper_get_layer, layer = x) |>
    purrr::reduce(cbind) |>
    `colnames<-`(tensor_coldata)

  final <- final[, !colSums(is.na(final)), drop = FALSE]

  names(final)
  return(final)
}





# sd ----------------------------------------------------------------------
message("Handling sd ...")

sd_coldata <- file_list_rds |>
  map(pluck, "image_file") |>
  unlist()

sd_bam_id <- file_list_rds |>
  map(pluck, "image_file") |>
  map(get_bam_id)

sd_rowdata <- file_list_rds |>
  map(pluck, "col_sd_df_final") |>
  map(pluck, "isize") |>
  pluck(1)

sd_assay <- file_list_rds |>
  map(pluck, "col_sd_df_final") |>
  map(pluck, "col_sd") |>
  map(as.matrix) |>
  setNames(sd_coldata) |>
  bind_cols() |>
  as.matrix() |>
  `dimnames<-`(list(sd_rowdata, sd_coldata))

sd_primary <- sd_bam_id |>
  map(lookup_primary,
    input_col = "bam_id",
    output_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()


sd_primary_boolean <- sd_bam_id |>
  map(lookup_primary_boolean,
    input_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()



sd_map <- tibble(
  primary = sd_primary,
  colname = sd_coldata[sd_primary_boolean]
)



# length ------------------------------------------------------------------
message("Handling length ...")

length_coldata <- file_list_rds |>
  map(pluck, "image_file") |>
  unlist()

length_bam_id <- file_list_rds |>
  map(pluck, "image_file") |>
  map(get_bam_id)

length_rowdata <- file_list_rds |>
  map(pluck, "col_sd_df_final") |>
  map(pluck, "isize") |>
  pluck(1)

length_assay <- file_list_rds |>
  map(pluck, "col_sd_df_final") |>
  map(pluck, "col_sum") |>
  map(as.matrix) |>
  setNames(length_coldata) |>
  bind_cols() |>
  as.matrix() |>
  `dimnames<-`(list(length_rowdata, length_coldata))



length_primary <- length_bam_id |>
  map(lookup_primary,
    input_col = "bam_id",
    output_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()


length_primary_boolean <- length_bam_id |>
  map(lookup_primary_boolean,
    input_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()

length_map <- tibble(
  primary = length_primary,
  colname = length_coldata[length_primary_boolean]
)


# sl ratio ----------------------------------------------------------------
message("Handling sl_ratio ...")

slratio_coldata <- file_list_rds |>
  map(pluck, "image_file") |>
  unlist()

slratio_bam_id <- file_list_rds |>
  map(pluck, "image_file") |>
  map(get_bam_id)


slratio_rowdata <- file_list_rds |>
  map(pluck, "bin_df_final") |>
  map(pluck, "bin") |>
  pluck(1)

slratio_assay <- file_list_rds |>
  map(pluck, "bin_df_final") |>
  map(pluck, "sl_ratio_corrected") |>
  map(as.matrix) |>
  setNames(slratio_coldata) |>
  bind_cols() |>
  as.matrix() |>
  `dimnames<-`(list(slratio_rowdata, slratio_coldata))

slratio_primary <- slratio_bam_id |>
  map(lookup_primary,
    input_col = "bam_id",
    output_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()

slratio_primary_boolean <- slratio_bam_id |>
  map(lookup_primary_boolean,
    input_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()


slratio_map <- tibble(
  primary = slratio_primary,
  colname = slratio_coldata[slratio_primary_boolean]
)


# cnv ---------------------------------------------------------------------
message("Handling cnv ...")

cnv_log2ratio_coldata <- file_list_rds |>
  map(pluck, "image_file") |>
  unlist()

cnv_log2ratio_bam_id <- file_list_rds |>
  map(pluck, "image_file") |>
  map(get_bam_id)


cnv_log2ratio_rowdata <- file_list_rds |>
  map(pluck, "bin_df_final") |>
  map(pluck, "bin") |>
  pluck(1)

cnv_log2ratio_assay <- file_list_rds |>
  map(pluck, "bin_df_final") |>
  map(pluck, "cnv_log2Ratio") |>
  map(as.matrix) |>
  setNames(cnv_log2ratio_coldata) |>
  bind_cols() |>
  as.matrix() |>
  `dimnames<-`(list(cnv_log2ratio_rowdata, cnv_log2ratio_coldata))



cnv_log2ratio_primary <- cnv_log2ratio_bam_id |>
  map(lookup_primary,
    input_col = "bam_id",
    output_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()

cnv_log2ratio_primary_boolean <- cnv_log2ratio_bam_id |>
  map(lookup_primary_boolean,
    input_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()


cnv_log2ratio_map <- tibble(
  primary = cnv_log2ratio_primary,
  colname = cnv_log2ratio_coldata[cnv_log2ratio_primary_boolean]
)

# motif -------------------------------------------------------------------
message("Handling motif_ratio ...")

motif_ratio_coldata <- file_list_rds |>
  map(pluck, "image_file") |>
  unlist()

motif_ratio_bam_id <- file_list_rds |>
  map(pluck, "image_file") |>
  map(get_bam_id)

motif_ratio_rowdata <- file_list_rds |>
  map(pluck, "bin_df_final") |>
  map(pluck, "bin") |>
  pluck(1)

motif_ratio_assay <- file_list_rds |>
  map(pluck, "bin_df_final") |>
  map(pluck, "s1_CT_ratio_corrected") |>
  map(as.matrix) |>
  setNames(motif_ratio_coldata) |>
  bind_cols() |>
  as.matrix() |>
  `dimnames<-`(list(motif_ratio_rowdata, motif_ratio_coldata))

motif_ratio_primary <- motif_ratio_bam_id |>
  map(lookup_primary,
    input_col = "bam_id",
    output_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()

motif_ratio_primary_boolean <- motif_ratio_bam_id |>
  map(lookup_primary_boolean,
    input_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()


motif_ratio_map <- tibble(
  primary = motif_ratio_primary,
  colname = motif_ratio_coldata[motif_ratio_primary_boolean]
)



# tensor ------------------------------------------------------------------
message("Handling tensor ...")
tensor_coldata <- file_list_rds |>
  map(pluck, "image_file") |>
  unlist()

tensor_bam_id <- file_list_rds |>
  map(pluck, "image_file") |>
  map(get_bam_id)


tensor_ratio_rowdata <- file_list_rds |>
  map(pluck, "tensor") |>
  map(dimnames) |>
  pluck(1)

tensor_assay <- file_list_rds |>
  map(pluck, "tensor") |>
  map(array_branch, margin = 3) |>
  lapply(FUN = map, .f = reshape2::melt) |>
  lapply(FUN = map, .f = tensor_helper)

# define a function to change NA to 0
replace_na_with_zero <- function(dataframe) {
  dataframe$value <- ifelse(is.na(dataframe$value), 0, dataframe$value)
  return(dataframe)
}

# change NA to 0
tensor_assay <- lapply(tensor_assay, function(sample) {
  lapply(sample, replace_na_with_zero)
})

layers <- tensor_assay[[1]] |>
  names()
names(layers) <- layers



tensor_assay_list <- map(layers,
  .f = tensor_helper_merge,
  tensor_assay = tensor_assay,
  tensor_coldata = tensor_coldata
)

# Define a function to replace NA with 0 in a matrix
replace_na_with_zero <- function(mat) {
  mat[is.na(mat)] <- 0L
  return(mat)
}

# Apply the function to each matrix in the list
# tensor_assay_list <- lapply(tensor_assay_list, replace_na_with_zero)


# tensor_assay_list <- bettermc::mclapply(X = layers,
#                         FUN = tensor_helper_merge,
#                         tensor_assay = tensor_assay,
#                         tensor_coldata = tensor_coldata,
#                         mc.cores = 2)


tensor_primary <- tensor_bam_id |>
  map(lookup_primary,
    input_col = "bam_id",
    output_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()

tensor_primary_boolean <- tensor_bam_id |>
  map(lookup_primary_boolean,
    input_col = "bam_id",
    meta_file = ulyses_meta
  ) |>
  unlist()


tensor_map <- tibble(
  primary = tensor_primary,
  colname = tensor_coldata[tensor_primary_boolean]
) |>
  list()


tensor_map_list <- rep(tensor_map, length(layers)) |>
  purrr::set_names(layers)


# mae ---------------------------------------------------------------------
message("Creating MAE object ...")
assay_list <- list(
  "sd" = sd_assay,
  "length" = length_assay,
  "sl_ratio" = slratio_assay,
  "cnv" = cnv_log2ratio_assay,
  "motif_ratio" = motif_ratio_assay
)

assay_list <- c(assay_list, tensor_assay_list)


map_list <- list(
  "sd" = sd_map,
  "length" = length_map,
  "sl_ratio" = slratio_map,
  "cnv" = cnv_log2ratio_map,
  "motif_ratio" = motif_ratio_map
)

map_list <- c(map_list, tensor_map_list) |>
  map(as.data.frame)


map_list_df <- listToMap(map_list)



ulyses_meta2 <- ulyses_meta |>
  column_to_rownames("bam_id")

mae <- MultiAssayExperiment(
  experiments = assay_list,
  colData = ulyses_meta2,
  sampleMap = map_list_df
)
message("MAE object saving ...")
saveRDS(object = mae, file = mae_filename)
message("MAE object saved to ", mae_filename, " ...")
