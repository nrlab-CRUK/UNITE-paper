
# ulyses main function
ulyses_main_gpu <- function(
    x,
    # motif parameters
    motif_type_vec = factor(c("u", "s", "umono", "smono"), 
      levels = c("u", "s", "umono", "smono")),
    motif_length_vec = 1:3,
    # functions
    bin_functions = "/home/nrlab/wang04/ulyses/explore_packaged/functions/bin_functions.r",
    frag_functions = "/home/nrlab/wang04/ulyses/explore_packaged/functions/frag_functions.r",
    overlap_functions = "/home/nrlab/wang04/ulyses/explore_packaged/functions/overlap_functions.r",
    
    # essential resource files 
    ref_csv_file = "/home/nrlab/wang04/ulyses/explore_packaged/resources/overlap_ref.csv",
    bins_anno_tiled_file = "/home/nrlab/wang04/ulyses/explore_packaged/resources/bins_anno_tiled.rds",
    
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
               "bin_mean_GC", 
               "bin_mean_mappability")
    
){
  
  ################################################################
  ## libs and sources
  ################################################################
  
  
  library(parallel, quietly = TRUE, warn.conflicts=FALSE)
  library(future, quietly = TRUE, warn.conflicts = FALSE)
  library(furrr)
  library(progressr)
  library(parallelly)
  #library(BiocGenerics)
  
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
  

# List of files to check
f_valid <- c(bin_functions, 
             frag_functions,
             overlap_functions,
             ref_csv_file,
             bins_anno_tiled_file
)

# Check if all files exist
if (!all(file.exists(f_valid))) {
  message("File(s) missing:")
  message(f_valid[!file.exists(f_valid)])
  stop("Error:Stopped due to missing file(s).")
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
      annotate_motif(genome = genome, 
        motif_type_vec = motif_type_vec,
        motif_length_vec = motif_length_vec ) %>%
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
      annotate_motif(genome = genome, 
        motif_length_vec = motif_length_vec,
        motif_type_vec = motif_type_vec)
    
    frag <- frag %>%
      plyranges::filter(frag_len >= fragQC_isize_min & frag_len <= fragQC_isize_max) %>%
      plyranges::filter(mapq >= 0) %>%
      plyranges::filter(frag_len >= isize_from & frag_len <= isize_to)
  }

  message("overlapping...")
  overlap_list <- parallel::mclapply(bin, plyranges::join_overlap_left, y = frag) %>%
    purrr::map(tibble::as_tibble) %>%
    purrr::imap(~ .x %>% dplyr::mutate(bin_size_name_tmp = .y))
    #my_table <- arrow::arrow_table(my_tibble)
return(overlap_list)

}
