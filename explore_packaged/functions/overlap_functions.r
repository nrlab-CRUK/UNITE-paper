

each_output_layer <- function(output_col, each_id_isize_combination_tibble, ref, ...){
  
  colname_to_use <- dplyr::filter(ref, output_col == {{output_col}}) %>% 
    dplyr::pull(source_col)
  
  
  method_to_use <- dplyr::filter(ref, output_col == {{output_col}}) %>% 
    dplyr::pull(method)
  
  #############################################################################
  # prepare functions
  #############################################################################
  
  motif_filter_nrow <- function(each_id_isize_combination_tibble, colname_to_use, output_col){
    
    # get motif
    motif <- stringr::str_split(output_col, "_")[[1]] %>% tail(n = 1)
    
    
    # get motif count result
    
    result <- each_id_isize_combination_tibble %>%
      dplyr::filter(.data[[colname_to_use]] %in% as.vector(motif)) %>%
      nrow() 
    
    #msg <- paste("Detected", result, motif, "...")
    #message(msg)
    return(result)
  }
  
  #############################################################################
  # calculate metrics based on method 
  #############################################################################
  
  
  result <- if(method_to_use == "sum") {
    sum(each_id_isize_combination_tibble %>% dplyr::pull({{colname_to_use}}), na.rm = TRUE)
  }else if(method_to_use == "mean"){
    mean(each_id_isize_combination_tibble %>% dplyr::pull({{colname_to_use}}), na.rm = TRUE)
  }else if(method_to_use == "motif_filter_nrow"){
    motif_filter_nrow(each_id_isize_combination_tibble, colname_to_use = colname_to_use, output_col = output_col)
  }else{
    stop("undefined calculation method, please check funcion `each_output_layer`")
  }
  
  return(result)
  
}





each_id_isize_combination <- function(x, output_layer, ref) {
  
  result <- sapply(output_layer, FUN = each_output_layer, each_id_isize_combination_tibble = x, ref = ref)
  
  return(result)
  
}

each_id_isize_combination_progressor <- function(x, output_layer, ref, p, ...){
  
  p()
  
  each_id_isize_combination(x = x, 
                            output_layer = output_layer,
                            ref = ref)
  
  
}


each_bin_size <- function(x, output_layer, ref, isize_from, isize_to, normalize_output = FALSE, ...){
  
  
  if(!tibble::is_tibble(x)){
    
    x <- tibble::as_tibble(x)
  }
  # report progress: which bin size 
  
  #bin_size_name <- unique(dplyr::pull(x, "bin_size_name_tmp"))
  #message(paste("Now processing...", bin_size_name))
  
  #############################################################################
  # compute output layers 
  #############################################################################
  
  id_isize_split <- x %>%
    dplyr::group_by(id, frag_len) %>%
    dplyr::group_split()
  
  
  keys <- x %>%
    dplyr::group_by(id, frag_len) %>%
    dplyr::group_keys()
  
  #p_each_bin_size <- progressr::progressor(steps = length(id_isize_split))
  
  result <- furrr::future_map(id_isize_split, 
                              each_id_isize_combination,
                              output_layer = output_layer, 
                              ref = ref
  ) %>%
    dplyr::bind_rows() %>%
    tibble::add_column(keys, .before = output_layer[1])
  
  
  #############################################################################
  # complete table (isize) and add missing values
  #############################################################################
  
  missing_frag_len <-setdiff( isize_from:isize_to, 
                              dplyr::filter(result, !is.na(frag_len)) %>% 
                                dplyr::pull(frag_len) %>% 
                                unique())
  
  
  fill_list <- ref %>%
    dplyr::filter(output_col %in% output_layer) %>%
    dplyr::select(output_col, null_filling_value)
  
  l <- as.list(dplyr::pull(fill_list, null_filling_value)) %>%
    setNames(dplyr::pull(fill_list, output_col))
  
  # crossing id and missing frag_len
  to_add <- tidyr::crossing(keys$id, missing_frag_len)
  
  names(to_add) <- c("id", "frag_len")
  
  # add the missing frag_len to id, fill values and drop the original frag_len 
  # with NA values
  
  result_tidy <- result %>%
    tibble::add_row(id = to_add$id, frag_len = to_add$frag_len ) %>%
    tidyr::complete(id, frag_len, fill = l) %>%
    tidyr::drop_na(frag_len) 
  
  # replace NAs for the GC and mappability. i.e. 'bin_status' layers
  
  if("bin_mean_GC" %in% names(result_tidy)) {
    
    gc_tibble <- result %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(bin_mean_GC = mean(bin_mean_GC))
    
    result_tidy <- result_tidy %>%
      dplyr::select(-bin_mean_GC)
    
    result_tidy <- dplyr::left_join(result_tidy, gc_tibble, by = "id")
    
  }
  
  
  if("bin_mean_mappability" %in% names(result_tidy)) {
    mappability_tibble <- result %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(bin_mean_mappability = mean(bin_mean_mappability))
    
    result_tidy <- result_tidy %>%
      dplyr::select(-bin_mean_mappability)
    
    result_tidy <- dplyr::left_join(result_tidy, mappability_tibble, by = "id")
  }
  
  #############################################################################
  # normalize layers
  #############################################################################
  
  min_max_norm <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  
  
  # generate norm groups
  
  layer_norm_groups <-  ref %>%
    dplyr::filter(output_col %in% output_layer) %>%
    dplyr::select(output_col, type_for_norm) %>%
    dplyr::group_by(type_for_norm) %>%
    dplyr::group_split(.keep = FALSE)
  
  # convert col to vector
  
  layer_norm_groups <- purrr::map(layer_norm_groups, .f = ~dplyr::pull(.))
  
  
  # create the must-have cols, i.e. id and frag_len
  
  id_len <- result_tidy %>% dplyr::select(id, frag_len)
  
  # norm all columns within the same 'type_for_norm' group together
  
  result_tidy_norm <- purrr::map(layer_norm_groups, 
                                 .f = ~ dplyr::select(result_tidy, all_of(.)) %>% min_max_norm()) %>%
    purrr::reduce(.f = tibble::tibble, .init = id_len)
  
  if (normalize_output) {
    return(result_tidy_norm)
  } else {
    return(result_tidy)
  }
}



each_bin_size_progressor <- function(x, output_layer, ref, isize_from, isize_to, p,  ...  ) {
  
  p()
  
  each_bin_size(x = x,
                output_layer = output_layer,
                ref = ref,
                isize_from = isize_from,
                isize_to = isize_to
  )
}





each_sample <- function(x, output_layer, ref, isize_from, isize_to, ...) {
  
  if(!is.list(x)){stop("Input must be a list.")}
  lapply(x, FUN = each_bin_size, output_layer = output_layer, ref = ref, isize_from = isize_from, isize_to = isize_to) %>%
    dplyr::bind_rows(.id = "tile_size")
  
}


# conver long tibble summary to image array
long_to_array <- function(obj) {
  
  channel_names <- unique(obj$channel)
  row_names <- unique(obj$id)
  col_names <- unique(obj$frag_len) %>% 
    paste("isize_", ., sep = "") 
  
  target_dim_names <- list(row_names, col_names, channel_names)
  
  target_dim <- c(length(row_names), 
                  length(col_names), 
                  length(channel_names))
  
  
  image_array <- base::array(obj$pixel, 
                             dim = target_dim, 
                             dimnames = target_dim_names)
  
  return(image_array)
}


# ulyses main function
ulyses_main <- function(
    x,
    
    # save or not
    save_olap = TRUE,
    save_summary_long = TRUE,
    # outputs
    bin_frag_overlap_file,
    summary_long_file,
    # functions
    bin_functions = "/home/nrlab/wang04/ulyses/explore_packaged/functions/bin_functions.r",
    frag_functions = "/home/nrlab/wang04/ulyses/explore_packaged/functions/frag_functions.r",
    overlap_functions = "/home/nrlab/wang04/ulyses/explore_packaged/functions/overlap_functions.r",
    
    # essential resource files 
    ref_csv_file = "/home/nrlab/wang04/ulyses/explore_packaged/resources/overlap_ref.csv",
    bins_anno_tiled_file = "/home/nrlab/wang04/ulyses/explore_QDNAseq_blacklist/bins_anno_tiled.rds",
    
    # frag annoation params
    which_genome = "hg19",
    fragQC_isize_min = 0,
    fragQC_isize_max = 2000,
    
    # bin-frag overlap params
    isize_from = 20,
    isize_to = 500,
    
    # image embedding params
    layers = c("n_isize",
               "n_motif_s1_C",
               "n_motif_s1_A",
               "n_motif_s1_G",
               "n_motif_s1_T",
               "n_neg_nm", 
               "n_pos_nm", 
               "bin_mean_GC", 
               "bin_mean_mappability")
    
){
  
  ################################################################
  ## libs and sources
  ################################################################
  
  
  library(parallel, quietly = TRUE, warn.conflicts=FALSE)
  library(bettermc)
  library(future, quietly = TRUE, warn.conflicts = FALSE)
  library(furrr)
  library(progressr)
  library(parallelly)
  
  library(tidyverse, quietly = TRUE, warn.conflicts=FALSE)
  library(Biobase, quietly = TRUE, warn.conflicts=FALSE)
  library(QDNAseq, quietly = TRUE, warn.conflicts=FALSE)
  library(Homo.sapiens, quietly = TRUE, warn.conflicts=FALSE)
  library(plyranges, quietly = TRUE, warn.conflicts=FALSE)
  library(GenomicRanges, quietly = TRUE, warn.conflicts=FALSE)
  library(GenomicAlignments, quietly = TRUE, warn.conflicts=FALSE)
  library(matrixStats, quietly = TRUE, warn.conflicts=FALSE)
  library(AnnotationHub, quietly = TRUE, warn.conflicts=FALSE)
  source(bin_functions)
  source(frag_functions)
  source(overlap_functions)
  
  ################################################################
  ## params QC
  ################################################################
  
  #### check if files exist
  f_valid <- c(bin_functions, 
               frag_functions,
               overlap_functions,
               ref_csv_file,
               bins_anno_tiled_file
  )
  
  if(!purrr::every(f_valid, file.exists)){
    message("File(s) missing:")
    message(f_valid[-which(file.exists(f_valid))])
    message("*******************************************")
    stop("Stopped due to missing file(s).")
  }
  
  #### read in resource files
  if (which_genome == "hg19") {
    library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
    genome <- BSgenome.Hsapiens.UCSC.hg19 
  } else if (which_genome == "hg38") {
    library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
    genome <- BSgenome.Hsapiens.UCSC.hg38 
  }
  
  ref <- readr::read_csv(ref_csv_file, show_col_types = FALSE)
  bin <- readRDS(bins_anno_tiled_file)
  
  if(length(bin) == 0 || is.na(bin)) {
    stop("bin annotation file lenght is 0, stop running!", call. = FALSE)
  }
  
  if(length(ref) == 0 || is.na(ref)) {
    stop("ref file lenght is 0, stop running!", call. = FALSE)
  }
  
  ################################################################
  ## function body 
  ################################################################
  
  if (!is(x, "GRanges")){
    message("Converting input to Fragment GRanges...")
    frag <- x %>%
      remove_outward_facing_readpairs() %>%
      curate_start_and_end() %>%
      `seqlengths<-`(seqlengths(genome)[levels(seqnames(.))]) %>%
      `genome<-`(seqinfo(genome)@genome[1])
    
    frag <- frag %>%
      remove_out_of_bound_reads() %>%
      annotate_motif(genome = genome) %>%
      plyranges::filter(frag_len >= fragQC_isize_min & frag_len <= fragQC_isize_max) %>%
      plyranges::filter(pos_mapq >= 0 & neg_mapq >= 0) %>%
      plyranges::filter(!stringr::str_detect(pos_cigar, pattern = "[ID]")) %>%
      plyranges::filter(!stringr::str_detect(neg_cigar, pattern = "[ID]")) %>%
      plyranges::filter(frag_len >= isize_from & frag_len <= isize_to)
    
  } else {
    
    frag <- x
    
    frag <- frag %>%
      remove_out_of_bound_reads() %>%
      plyranges::mutate(frag_len = BiocGenerics::width(frag), frag_len_seen = 1) %>%
      annotate_motif(genome = genome, motif_length_vec = c(1, 2, 3, 4))
    
    frag <- frag %>%
      plyranges::filter(frag_len >= fragQC_isize_min & frag_len <= fragQC_isize_max) %>%
      plyranges::filter(mapq >= 0) %>%
      plyranges::filter(frag_len >= isize_from & frag_len <= isize_to)
  }
  
  
  message("overlapping...")
  overlap_list <- bettermc::mclapply(bin, plyranges::join_overlap_left, y = frag) %>%
    purrr::map(tibble::as_tibble) %>%
    purrr::imap(~ .x %>% dplyr::mutate(bin_size_name_tmp = .y))
  
  # save file or not
  if (save_olap) {
    saveRDS(overlap_list, bin_frag_overlap_file)
  }
  
  # set up parallelization 
  
  plan(list(
    tweak(multisession, workers = parallelly::availableCores(method = "system") %/% 32),
    tweak(multisession, workers = 32)))
  
  message("summarizing metrics...")
  summary_list <- furrr::future_map(overlap_list, 
                                    each_bin_size, 
                                    output_layer = layers,
                                    ref = ref,
                                    isize_from = isize_from,
                                    isize_to = isize_to
  )
  
  message("summarizing long table...")
  tibble_long <- export_summary_table_long(x = summary_list, 
                                           layers = layers)
  
  # calculate ulyses object
  tibble_long_grouped <- tibble_long %>% 
    # the order is important
    dplyr::arrange(tile_size, channel,frag_len, id) %>%
    dplyr::group_by(tile_size)
  
  
  # save file or not
  if (save_summary_long) {
    saveRDS(tibble_long_grouped, summary_long_file)
  }
  
  tibble_long_grouped_keys <- dplyr::group_keys(tibble_long_grouped)
  
  tibble_long_grouped_list <- tibble_long_grouped %>% 
    dplyr::group_split() %>%
    magrittr::set_names(tibble_long_grouped_keys$tile_size)
  
  message("converting long tibble to images...")
  image_list <- purrr::map(tibble_long_grouped_list, long_to_array)
  
  message("Done...")
  return(image_list)
  
}

