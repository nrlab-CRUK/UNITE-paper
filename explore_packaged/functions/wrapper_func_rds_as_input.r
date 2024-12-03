

ulyses_wrapper <- function(rdsfile, 
	block_size = NA,
	slurm_partition = NA,
	package_path = NA
) {


if(is.null(slurm_partition) || is.na(slurm_partition)) {
  slurm_partition = "epyc"
}

if(is.null(package_path) || is.na(package_path)) {
  package_path = "/home/nrlab/wang04/ulyses/explore_packaged" 
}


# output file
output_path = dirname(rdsfile)
output_prefix = basename(rdsfile)
image_file = file.path(output_path, paste(output_prefix, ".ulyses_image.rds", sep = ''))
bin_frag_overlap_file = file.path(output_path, paste(output_prefix, ".bin_frag_olap.rds", sep = ''))
summary_long_file = file.path(output_path, paste(output_prefix, ".ulyses_summary_long.rds", sep = '')) 

save_olap = FALSE
save_summary_long = FALSE

# rslurm params
rslurm_jobname = paste(output_prefix, "_ulyses", sep = "")
rslurm_options = list(time = '10:00:00',  mem = '128G' , partition = slurm_partition)


# functions
bin_functions = file.path(package_path, "functions/bin_functions.r")
frag_functions = file.path(package_path, "functions/frag_functions.r")
overlap_functions = file.path(package_path, "functions/overlap_functions.r")

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

# image embedding params
layers = c("n_isize",
           "n_motif_s1_C",
           "n_motif_s1_A",
           "n_motif_s1_G",
           "n_motif_s1_T",
           "n_motif_s4_CCCA",
           "bin_mean_GC", 
           "bin_mean_mappability")

###############################################################################
# libs 
###############################################################################

library(tidyverse, quietly = TRUE, warn.conflicts=FALSE)
library(Rsamtools, quietly = TRUE, warn.conflicts=FALSE)
library(rslurm, quietly = TRUE, warn.conflicts=FALSE)
library(abind)
library(progressr)
source(bin_functions)
source(frag_functions)
source(overlap_functions)

###############################################################################
# preprossing
###############################################################################


#rslurm tmp folder
mainDir <- dirname(rdsfile) 
subDir <-  paste0(basename(rdsfile), "_TEMP_DIR_debug")
tmpDir <- file.path(mainDir, subDir) 

dir.create(tmpDir, showWarnings = FALSE)
setwd(tmpDir)


ref <- readr::read_csv(ref_csv_file, show_col_types = FALSE)
bin <- readRDS(bins_anno_tiled_file)

if(length(bin) == 0 || is.na(bin)) {
			stop("bin annotation file lenght is 0, stop running!", call. = FALSE)
}

if(length(ref) == 0 || is.na(ref)) {
			stop("ref file lenght is 0, stop running!", call. = FALSE)
}


###############################################################################
# Splitting rds into chunks 
###############################################################################

message(paste("Splitting GRanges into chunks, each chunk has ", block_size, " fragments."))


  frag <- readRDS(rdsfile)

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


###############################################################################
# utilize rslurm to parallelize the summary process
###############################################################################
message("Submitting overlap summary jobs via rslurm.")

job <- rslurm::slurm_map(bam_chunk_list,
                         f = ulyses_main,
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
                         layers = layers, 
                         save_olap = save_olap,
                         save_summary_long = save_summary_long,
                         bin_frag_overlap_file = bin_frag_overlap_file,
                         summary_long_file = summary_long_file,
                         jobname = rslurm_jobname,
                         nodes = 1000, 
                         cpus_per_node = 1,
                         processes_per_node = 1,
                         slurm_options = rslurm_options,
                         submit = TRUE)

# delete the bam_chunk_list to save MEM
rm(bam_chunk_list)

message("Waiting for rslurm output files... \n")
image_chunks <- get_slurm_out(job, outtype = 'raw', wait = TRUE)
message("Rslurm output files gathered! \n")
setwd(mainDir)
#rslurm::cleanup_files(job)

# merge into single image object 


message("Gathered slurm output, now calculate combined image array...")
f1 <- function(x, bin_size, layer){
	x[[bin_size]][,,layer]
	
}

f2 <- function(bin_size, file, layer){

	lapply(file, f1, bin_size = bin_size, layer = layer )

	

}

f3 <- function(file, bin_size, layer) {

  lapply(bin_size, f2, file, layer = layer)
}

matrix_reduce_mean <- function(x) {Reduce(`+`, x) / length(x)}
matrix_reduce_sum <- function(x) {Reduce(`+`, x)}

tmp <- purrr::map(as.list(layers) %>% setNames(layers), 
                  .f = f3, 
                  bin_size = 1:length(bin), 
                  file = image_chunks)


newlist <- list()
for (i in names(tmp)) {

   i_method <- dplyr::filter(ref, output_col == !!i) %>% dplyr::pull(method) %>% unique()


   if (i_method == "mean") {
      r <- purrr::map(tmp[[i]], .f = matrix_reduce_mean) %>% list()

   } else if(i_method %in% c("sum", "motif_filter_nrow")) {
      r <- purrr::map(tmp[[i]], .f = Reduce, f = `+`) %>% list()
   }

  newlist <- append(newlist, r)


}

newlist2 <- list()
for (bin_id in 1:length(bin)){


  r <- purrr::map(newlist, bin_id) %>% list()
  newlist2 <- append(newlist2, r)
}


final_image <- purrr::map(newlist2, abind::abind, along = 3) %>%
  setNames(names(bin))

for (bin_name in names(bin)) {

  dimnames(final_image[[bin_name]]) <- dimnames(image_chunks[[1]][[bin_name]])

}

message("Saving image file:")
message(image_file)
saveRDS(object = final_image, image_file )

#message("Deleting name sorted bam file")
#base::unlink(name_sort_bam_filename)

message("Done")
}
