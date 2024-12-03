


ulyses_wrapper <- function(bamfile, 
                           bamfile_index = NA,
                           block_size = NA,
                           slurm_partition = NA,
                           package_path = NA,
                           input_file_type = NA
) {
  
  
  if(is.null(slurm_partition) || is.na(slurm_partition)) {
    
    slurm_partition = "epyc"
    
  }
  
  if(is.null(package_path) || is.na(package_path)) {
    package_path = "/home/nrlab/wang04/ulyses/explore_packaged" 
    
  }
  

  if(is.null(input_file_type) || is.na(input_file_type)) {
    # detect suffix of bamfile
    bamfile_suffix <- tools::file_ext(bamfile)
    if(bamfile_suffix %in% c("bam", "BAM", "Bam")) {
      input_file_type = "bam"
    } else if(bamfile_suffix %in% c("rds", "RDS", "Rds")) {
      input_file_type = "rds"
    } else {
      stop("input_file_type is not specified, please check the input bamfile or rdsfile.")
    }
    
  }
  
  if(is.na(bamfile_index) && input_file_type == "bam") {
    bamfile_index = bamfile 
    
  }
  
  # output file
  output_path = dirname(bamfile)
  output_prefix = basename(bamfile)
  # a bin_size_list containing id frag_len overlaped 
  overlap_file = file.path(output_path, paste(output_prefix, ".ulyses_bind_olap_chunks.rds", sep = ''))
  # a bin_size_list containing layer/channel  summary
  layer_summary_file = file.path(output_path, paste(output_prefix, ".ulyses_layer_summary.rds", sep = ''))
  # a bin size list containing tensors
  image_file = file.path(output_path, paste(output_prefix, ".ulyses_image.rds", sep = ''))
  
  
  # rslurm params
  rslurm_jobname = paste(output_prefix, "_ulyses", sep = "")
  rslurm_options = list(time = '10:00:00',  mem = '16G' , partition = slurm_partition)
  
  
  # functions
  bin_functions = file.path(package_path, "functions/bin_functions.r")
  frag_functions = file.path(package_path, "functions/frag_functions.r")
  overlap_functions = file.path(package_path, "functions/overlap_functions.r")
  ulyses_main_gpu_functions = file.path(package_path, "functions/ulyses_main_gpu.r")
  
  # essential resource files 
  ref_csv_file = file.path(package_path, "resources/overlap_ref.csv")
  bins_anno_tiled_file = file.path(package_path, "resources/bins_anno_tiled.rds")
  
  # frag annoation params
  which_genome = "hg19"
  fragQC_isize_min = 0
  fragQC_isize_max = 2000
  
  # bin-frag overlap params
  isize_from = 20
  isize_to = 500

  # set up motif params: this will affect what cols in the overlap table
    motif_type_vec = factor(c("u", "s", "umono", "smono"), levels = c("u", "s", "umono", "smono"))

    motif_length_vec = 1:3

  
  # set up layers 

  if(input_file_type == "rds"){
    layers = c(
           "n_isize",
            "n_motif_umono3_A",
            "n_motif_umono3_C",
            "n_motif_umono3_G",
            "n_motif_umono3_T",
            "n_motif_umono2_A",
            "n_motif_umono2_C",
            "n_motif_umono2_G",
            "n_motif_umono2_T",
           "n_motif_umono1_A",
            "n_motif_umono1_C",
            "n_motif_umono1_G",
            "n_motif_umono1_T",
            "n_motif_smono1_A",
            "n_motif_smono1_C",
            "n_motif_smono1_G",
            "n_motif_smono1_T",
            "n_motif_smono2_A",
            "n_motif_smono2_C",
            "n_motif_smono2_G",
            "n_motif_smono2_T",
            "n_motif_smono3_A",
            "n_motif_smono3_C",
            "n_motif_smono3_G",
            "n_motif_smono3_T",
           "n_motif_s1_C",
           "n_motif_s1_A",
           "n_motif_s1_G",
           "n_motif_s1_T",
           "bin_mean_GC", 
           "bin_mean_mappability")
  } else if(input_file_type == "bam"){
    layers = c(
           "n_isize",
            "n_motif_umono3_A",
            "n_motif_umono3_C",
            "n_motif_umono3_G",
            "n_motif_umono3_T",
            "n_motif_umono2_A",
            "n_motif_umono2_C",
            "n_motif_umono2_G",
            "n_motif_umono2_T",
           "n_motif_umono1_A",
            "n_motif_umono1_C",
            "n_motif_umono1_G",
            "n_motif_umono1_T",
            "n_motif_smono1_A",
            "n_motif_smono1_C",
            "n_motif_smono1_G",
            "n_motif_smono1_T",
            "n_motif_smono2_A",
            "n_motif_smono2_C",
            "n_motif_smono2_G",
            "n_motif_smono2_T",
            "n_motif_smono3_A",
            "n_motif_smono3_C",
            "n_motif_smono3_G",
            "n_motif_smono3_T",
             "n_motif_s1_C",
             "n_motif_s1_A",
             "n_motif_s1_G",
             "n_motif_s1_T",
             "n_neg_nm", 
             "n_pos_nm", 
             "bin_mean_GC", 
             "bin_mean_mappability")
  }

  # extract the motif_type from layers using regex, n_motif_umonon3_A -> umono, n_motif_s10_AAAAAAAAAA -> s
  motif_type_vec_used <- layers %>% 
    purrr::keep(~grepl("n_motif_", .x)) %>%
    stringr::str_extract(pattern = "(?<=n_motif_)(.*?)(?=\\d)") %>% 
    unique()
  # extract the motif_length from layers using regex.
  motif_length_vec_used <- layers %>% 
    purrr::keep(~grepl("n_motif_", .x)) %>%
    stringr::str_extract(pattern = "\\d+") %>% 
    as.integer() %>% 
    unique()
  # stop if no all of the motif_type_vec_used are in motif_type_vec
  if(!all(motif_type_vec_used %in% motif_type_vec)) {
    stop("motif_type_vec_used is not a subset of motif_type_vec, please check layers and motif_type_vec.")
  }
  # stop if no all of the motif_length_vec_used are in motif_length_vec
  if(!all(motif_length_vec_used %in% motif_length_vec)) {
    stop("motif_length_vec_used is not a subset of motif_length_vec, please check layers and motif_length_vec.")
  }

  
  ###############################################################################
  # libs 
  ###############################################################################
  
  library(tidyverse, quietly = TRUE, warn.conflicts=FALSE)
  library(Rsamtools, quietly = TRUE, warn.conflicts=FALSE)
  library(rslurm, quietly = TRUE, warn.conflicts=FALSE)
  library(abind)
  library(tictoc)
  library(progressr)
  library(patchwork)
  library(EBImage)
  library(gridExtra)
  library(grid)
  library(reshape2)
  library(purrr,  warn.conflicts = FALSE)
  library(parallel,  warn.conflicts = FALSE)
  library(fastDummies)

  source(bin_functions)
  source(frag_functions)
  source(overlap_functions)
  source(ulyses_main_gpu_functions)
  

###############################################################################
# functions
###############################################################################

motif_filter_nrow_rebuild <- function(df, colname_to_use, output_col) {
  # Only keep columns: id, frag_len, colname_to_use
  #df <- df[, c('id', 'frag_len', colname_to_use), drop = FALSE]
  
  df <- df %>% dplyr::select(id, frag_len, !!colname_to_use)
  
  # add prefix 'n_' to colname_to_use
  new_colname_to_use <- paste0('n_', colname_to_use)
  
  # rename colname_to_use to new_colname_to_use
  colnames(df)[colnames(df) == colname_to_use] <- new_colname_to_use
  
  # Create dummy variables using model.matrix
  df_dummy <- fastDummies::dummy_cols(df, 
  select_columns = new_colname_to_use, 
  remove_selected_columns = TRUE, 
  ignore_na = TRUE)
  
  
  result <- df_dummy %>%
    dplyr::group_by(id, frag_len, .drop = FALSE) %>%
    dplyr::summarise(!!output_col := sum(.data[[output_col]]), .groups = 'drop') %>%
    arrange(id, frag_len)
  
  return(result)
}




each_output_layer_rebuild <- function(output_col, input_df, ref, ...){
  
  # report output_col
  message(paste("Handling layer:", output_col, "\n"))
  
  # Filter reference data frame based on output_col
  colname_to_use <- ref %>%
    filter(output_col == !!output_col) %>%
    pull(source_col)
  
  # if colname_to_use doesn't exist in input_df, return NULL
  if (!colname_to_use %in% colnames(input_df)) {
    message(paste("colname_to_use:", colname_to_use, "doesn't exist in input_df, return NULL."))
    result <- NULL
  } else {
  
  method_to_use <- ref %>%
    filter(output_col == !!output_col) %>%
    pull(method) 
  
  # Calculate metrics based on method
  if (method_to_use %in% c('sum')) {
    result <- input_df %>%
      select(id, frag_len, !!colname_to_use) %>%
      as_tibble() %>%
      group_by(id, frag_len, .drop = FALSE) %>%
      summarise(!!output_col := sum(.data[[colname_to_use]], na.rm = TRUE), .groups = "drop") %>%
      arrange(id, frag_len)
  } else if (method_to_use == 'mean') {
    result <- input_df %>%
      select(id, frag_len, !!colname_to_use) %>%
      group_by(id, frag_len, .drop = FALSE) %>%
      summarise(!!output_col := mean(.data[[colname_to_use]], na.rm = TRUE), .groups = "drop") %>%
      arrange(id, frag_len)
  } else if (method_to_use == 'motif_filter_nrow') {
    result <- motif_filter_nrow_rebuild(df = input_df, colname_to_use = colname_to_use, output_col = output_col)
  } else {
    stop("Undefined calculation method, please check function `each_output_layer_rebuild` and `ref` data frame.")
  }
  }
  
  return(result)
}


each_bin_size_processing <- function(bin_size, combined, layers, ref, isize_from, isize_to) {
  cat("Handling bin_size:", bin_size, "\n")
  
  df1 <- combined[[bin_size]] %>%
    mutate(
      id = factor(id, levels = unique(id), ordered = TRUE),
      frag_len = factor(frag_len, levels = seq(isize_from, isize_to), ordered = TRUE)
    )
  
  # Report nlevels of id and frag_len
  message(paste("nlevels of id:", nlevels(df1$id), "\n"))
  message(paste("nlevels of frag_len:", nlevels(df1$frag_len), "\n"))
  
  # Map each output_col to the result of each_output_layer_rebuild function
  layers_result <- lapply(layers, each_output_layer_rebuild, input_df = df1, ref = ref)
  
  # Merge all DataFrames in the list using an outer join
  merged_df <- purrr::reduce(layers_result, ~dplyr::full_join(.x, .y, by = c('id', 'frag_len'))) %>%
    # Convert id and frag_len to factor, without order
    mutate(
      id = factor(id, levels = unique(id), ordered = FALSE),
      frag_len = factor(frag_len, levels = seq(isize_from, isize_to), ordered = FALSE)
    )
  
  return(merged_df)
}




  ###############################################################################
  # preprossing
  ###############################################################################
  
  
  #rslurm tmp folder
  mainDir <- dirname(bamfile) 
  subDir <-  paste0(basename(bamfile), "_TEMP_DIR_debug")
  tmpDir <- file.path(mainDir, subDir) 
  
  dir.create(tmpDir, showWarnings = FALSE)
  setwd(tmpDir)
  
  
  ref <- readr::read_csv(ref_csv_file, show_col_types = FALSE)
  bin <- readRDS(bins_anno_tiled_file)
  
  
  ###############################################################################
  # Splitting bam into chunks 
  ###############################################################################
  
  #name sorted bam file name
  #name_sort_bam_filename <- paste(bamfile, ".name_sorted.bam", sep = "")
  # name sort bam if it doesn't exist
  #if(!file.exists(name_sort_bam_filename)) {
  #  Rsamtools::sortBam(bamfile, byQname = TRUE, maxMemory = 1024, destination =  paste(bamfile, ".name_sorted", sep = "")) 
  #}
  
  # if input_file_type is bam , excute below
if(input_file_type == "bam") {
  message(paste("Splitting bam file into chunks, size = ", block_size, "fragments."))
  
  bam_fl <- Rsamtools::BamFile(bamfile, index = bamfile_index, asMates = TRUE)
  open(bam_fl)
  #bam_count <- countBam(bam_fl)$records
  yieldSize(bam_fl) <- block_size
  bam_chunk_list <- list()
  
  repeat {
    obj <- bam_to_galp(bam_fl) 
    obj_len <- length(obj)
    obj_l <- list(obj)
    if (obj_len == 0) {close(bam_fl) ; break}
    bam_chunk_list <- append(bam_chunk_list, obj_l) 
    message(paste(obj_len, "fragments retrieved..."))
  }
} else if(input_file_type == "rds") {
  message(paste("Splitting rds file into chunks, size = ", block_size, "fragments."))
  frag <- readRDS(bamfile)

boundaries <- seq(from = 0, to = length(frag), by = block_size) %>% 
  as.integer()

if(boundaries[length(boundaries)] < length(frag)) {
  boundaries <- append(boundaries, length(frag))
}

for (i in seq(1, length(boundaries))) {

  if(i == 1){
    bam_chunk_list <- list(frag[1:boundaries[i+1]])
  } else if (i < length(boundaries) & i > 1) {
    index <- seq(from = boundaries[i] +1 , to = boundaries[i + 1], by = 1)
    obj <- frag[index]
    bam_chunk_list <- append(bam_chunk_list, obj)

  }


}
}

  
  ###############################################################################
  # utilize rslurm to parallelize the summary process
  ###############################################################################
  
  message("Submitting overlap summary jobs via rslurm.")
  
  job <- rslurm::slurm_map(bam_chunk_list,
                           f = ulyses_main_gpu,
                           bin_functions = bin_functions,
                           frag_functions = frag_functions,
                           overlap_functions = overlap_functions,
                           ref_csv_file = ref_csv_file,
                           bins_anno_tiled_file = bins_anno_tiled_file,
                           which_genome = which_genome,
                           fragQC_isize_min = fragQC_isize_min,
                           fragQC_isize_max = fragQC_isize_max,
                           isize_from = isize_from,
                           isize_to = isize_to,
                           motif_type_vec = motif_type_vec,
                           motif_length_vec = motif_length_vec,
                           layers = layers, 
                           jobname = rslurm_jobname,
                           nodes = 1000, 
                           cpus_per_node = 1,
                           processes_per_node = 1,
                           slurm_options = rslurm_options,
                           #global_objects = global_objects,
                           submit = TRUE)
  
  # delete the bam_chunk_list to save MEM
  rm(bam_chunk_list)
  
  message("Waiting for rslurm output files... \n")
  overlap_chunks <- get_slurm_out(job, outtype = 'raw', wait = TRUE)
  message("Rslurm output files gathered! \n")
  setwd(mainDir)
  #rslurm::cleanup_files(job)
  
  # merge into single image object 
  # bind rows in each list element
  combined <- overlap_chunks %>% purrr::list_transpose() %>% purrr::map(bind_rows)

  # remove any rows where id or frag_len is NA, split into two steps to avoid error
  combined <- combined %>% purrr::map(~filter(., !is.na(id)))
  combined <- combined %>% purrr::map(~filter(., !is.na(frag_len)))
  
  message("Saving olap file:")
  message(overlap_file)
  saveRDS(object = combined, overlap_file)
  
  ###############################################################################
  # summarise overlap table
  ###############################################################################
  # Get unique bin_size levels
  bin_size_levels <- names(combined)
  
  # Initialize the summary_list list
  summary_list <- parallel::mclapply(bin_size_levels, 
                           each_bin_size_processing, 
                           combined = combined, 
                           layers = layers, 
                           ref = ref, 
                           isize_from = isize_from, 
                           isize_to = isize_to, 
                           mc.cores = parallel::detectCores())
  
  # set the names of the list
  names(summary_list) <- bin_size_levels 
  
  message("Saving summary_list:", layer_summary_file)
  saveRDS(object = summary_list, layer_summary_file)
  
  
  ###############################################################################
  # convert summary table long and then to ulyses object (tensors) 
  ###############################################################################
  
  tibble_long <- export_summary_table_long(x = summary_list, layers = layers)
  
  
  # calculate ulyses object
  tibble_long_grouped <- tibble_long %>% 
    # the order is important
    dplyr::arrange(tile_size, channel,frag_len, id) %>%
    dplyr::group_by(tile_size)
  tibble_long_grouped_keys <- dplyr::group_keys(tibble_long_grouped)
  tibble_long_grouped_list <- tibble_long_grouped %>% 
    dplyr::group_split() %>%
    magrittr::set_names(tibble_long_grouped_keys$tile_size)
  message("converting long tibble to images...")
  image_list <- purrr::map(tibble_long_grouped_list, long_to_array)


  message("Saving image file:")
  message(image_file)
  saveRDS(object = image_list, image_file )
  message("Saved.")

  ##############################################################################
  # plot ulyses images
  ##############################################################################

}
