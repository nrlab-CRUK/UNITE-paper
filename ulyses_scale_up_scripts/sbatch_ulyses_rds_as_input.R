library(rslurm)
library(tidyverse)

###############################################################################
# params
###############################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least rdsfile_path must be supplied (input file). \n", call.=FALSE)
}

# rdsfile
rdsfile_path <-  as.character(args[1])
slurm_partition <- as.character(args[2])
package_path <- as.character(args[3])
block_size <- as.numeric(args[4]) 


if(is.null(rdsfile_path) || is.na(rdsfile_path)) {
  rdsfile_path <- "." 
  
}

if(is.null(slurm_partition) || is.na(slurm_partition)) {
	if(stringr::str_detect(rdsfile_path, pattern = "scratchc")){
  		slurm_partition <-  as.character("epyc")
	}else if(stringr::str_detect(rdsfile_path, pattern = "scratcha")){
  		slurm_partition <-  as.character("general")
	}else if(stringr::str_detect(rdsfile_path, pattern = "scratchb")){
  		slurm_partition <-  as.character("general")
	} else {
		stop("Cannot infer which partition from your rdsfile_path, please indicate in parameter settings.")
	}
  
}

if(is.null(package_path) || is.na(package_path)) {
  package_path <- as.character("/home/nrlab/wang04/ulyses/explore_packaged")
  
}

if(is.null(block_size) || is.na(block_size)) {
	block_size <- as.numeric(1000000)
  
}

# get the files

rdsfile <- list.files(path = rdsfile_path, 
	pattern = "*.rds$", 
	full.names = TRUE)


run_file_list <- tibble(rdsfile = rdsfile, 
	block_size = block_size,
	slurm_partition = slurm_partition,
	package_path = package_path)

# filter by existence of rds files

run_file_list_filtered <- dplyr::filter(run_file_list, 
	!file.exists(paste(rdsfile, ".ulyses_image.rds", sep = "")))

# rslurm function

source("/home/nrlab/wang04/ulyses/explore_packaged/functions/wrapper_func_rds_as_input.r")

# run rslurm_apply on the filtered list

sjob <- rslurm::slurm_apply(
	f = ulyses_wrapper, 
	params = run_file_list_filtered,
	pkgs = c("tidyverse", "rslurm"),
	rscript_path = "/home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript",
	job_array_task_limit = 10,
	nodes = 3000, 
        cpus_per_node = 1,
        processes_per_node = 1,
	jobname = "u-hq-30d",
	slurm_options = list(time = '30-0',  mem = '1G' , partition = slurm_partition),
	submit = TRUE)


message("Mission Completed, Goodbye!")