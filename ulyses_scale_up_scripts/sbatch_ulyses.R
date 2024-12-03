library(rslurm)
library(tidyverse)
library(argparser)
library(stringr)

# Define command-line arguments
parser <- arg_parser(description = "Run ulyese for all files in a folder.")

parser <- add_argument(parser, "--bamfile_path", type = "character",
                       help = "Path to the BAM file.", default = ".")
parser <- add_argument(parser, "--slurm_partition", type = "character",
                       help = "Slurm partition.", default = "infer")
parser <- add_argument(parser, "--package_path", type = "character",
                       help = "Package path.",
                       default = "/home/nrlab/wang04/ulyses/explore_packaged")
parser <- add_argument(parser, "--block_size", type = "numeric",
                       help = "Block size.", default = 2000000)
parser <- add_argument(parser, "--input_file_type", type = "character",
                       help = "Input file type.", default = "bam")
parser <- add_argument(parser, "--ulyses_wrapper", type = "character",
                       help = "Ulyses wrapper function.",
                       default = "/home/nrlab/wang04/ulyses/explore_packaged/functions/wrapper_func4.R")
parser <- add_argument(parser, "--slurm_time", type = "character",
                       help = "Slurm time.", default = "0-5")
parser <- add_argument(parser, "--slurm_mem", type = "character",
                       help = "Slurm mem.", default = "32G")
parser <- add_argument(parser, "--clear_cache", type = "boolean",
                       help = "Clear slurm cache.", default = FALSE)
parser <- add_argument(parser, "--slurm_job_array_task_limit", type = "numeric",
                       help = "Run x samples in parallel.", default = 50)

# add parameter, replace output or not, default is FALSE.
parser <- add_argument(parser, "--replace_output", type = "boolean",
                       help = "Replace output or not.", default = FALSE)




# Parse command-line arguments
args <- parse_args(parser)

# Extract values from parsed arguments
bamfile_path <- args$bamfile_path
slurm_partition <- args$slurm_partition
package_path <- args$package_path
block_size <- args$block_size
input_file_type <- args$input_file_type
ulyses_wrapper <- args$ulyses_wrapper
slurm_time <- args$slurm_time
slurm_mem <- args$slurm_mem
clear_cache <- args$clear_cache
slurm_job_array_task_limit <- args$slurm_job_array_task_limit
replace_output <- args$replace_output


if(slurm_partition == "infer") {
	if(stringr::str_detect(bamfile_path, pattern = "scratchc")){
  		slurm_partition <-  as.character("epyc")
	}else if(stringr::str_detect(bamfile_path, pattern = "scratcha")){
  		slurm_partition <-  as.character("general")
	}else if(stringr::str_detect(bamfile_path, pattern = "scratchb")){
  		slurm_partition <-  as.character("general")
	} else {
		stop("Cannot infer which partition from your bamfile_path, please indicate in parameter settings.")
	}
  
}



# get the files

if(input_file_type == "bam"){

	bamfile <- list.files(path = bamfile_path, 
	pattern = "*.bam$", 
	full.names = TRUE)

	bamfile_ids <- gsub(".bam$", "", bamfile)


	bamfile_index <- list.files(path = bamfile_path, 
	pattern = "*.bai$",
	full.names = TRUE)

	# make bamfile_ids, if detected ".bam.bai" in the bamfile_index, remove it.
	# else if detected ".bai" in the bamfile_index, remove it.
	# else stop.
	if(length(grep(".bam.bai$", bamfile_index)) > 0){
		bamfile_index_ids <- gsub(".bam.bai$", "", bamfile_ids)
	} else if(length(grep(".bai$", bamfile_index)) > 0){
		bamfile_index_ids <- gsub(".bai$", "", bamfile_ids)
	} else {
		stop("Cannot find bamfile_index, please check your bamfile_path.")
	}

	# if any unique(bamfile_ids) not in bamfile_index_ids, stop.
	if(length(setdiff(unique(bamfile_ids), unique(bamfile_index_ids))) > 0){
		# report which bamfile_ids not in bamfile_index
		message("The following bamfile_ids not in bamfile_index: ",
		paste(setdiff(unique(bamfile_ids), unique(bamfile_index_ids)), collapse = ", "))
		# stop
		stop("Cannot find bamfile_index, please check your bamfile_path.")
	}


} else if(input_file_type == "rds"){
	
	bamfile <- list.files(path = bamfile_path, 
	pattern = "*.frag.tsv.bgz.GRanges.rds.1M.rds$", 
	full.names = TRUE)

	bamfile_index <- list.files(path = bamfile_path, 
	pattern = "*.frag.tsv.bgz.GRanges.rds.1M.rds$",
	full.names = TRUE)

}



bamfile_index <- gsub(".bai$", "", bamfile_index)

run_file_list <- tibble(bamfile = bamfile, 
	bamfile_index = bamfile_index,
	block_size = block_size,
	slurm_partition = slurm_partition,
	package_path = package_path)

# filter by existence of rds files

if(!replace_output){
run_file_list_filtered <- dplyr::filter(run_file_list, 
	!file.exists(paste(bamfile, ".ulyses_image.rds", sep = "")))

} else {
run_file_list_filtered <- run_file_list
}

# check if there is any file to run
if(nrow(run_file_list_filtered) == 0){
	message("No file to run, Goodbye!")
	quit(save = "no", status = 0)
}

# report how many files to run
message("There are ", nrow(run_file_list_filtered), " files to run.")
# report the number of files filtered
message("There are ", nrow(run_file_list) - nrow(run_file_list_filtered), " files won't run because their ulyses_image.rds files exist.")

# the ulyses_wrapper function sourced:
source(ulyses_wrapper)

sjob <- rslurm::slurm_apply(
	f = ulyses_wrapper, 
	params = run_file_list_filtered,
	pkgs = c("tidyverse", "rslurm"),
	rscript_path = "/home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript",
	job_array_task_limit = slurm_job_array_task_limit,
	nodes = 3000, 
        cpus_per_node = 1,
        processes_per_node = 1,
	jobname = "ulyses",
	slurm_options = list(time = slurm_time,  mem = slurm_mem , partition = slurm_partition),
	submit = TRUE)

message("Mission Completed, Goodbye!")
# clean up if clear_cache is TRUE
if(clear_cache == "TRUE"){
	cleanup_files(sjob)
}

