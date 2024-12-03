library(rslurm)
library(tidyverse)

###############################################################################
# params
###############################################################################
args <- commandArgs(trailingOnly=TRUE)
# stop if no args
if (length(args)==0) {
  stop("At least bamfile_path must be supplied (input file). \n", call.=FALSE)
}

# bamfile
bamfile_path <-  as.character(args[1])
slurm_partition <- as.character(args[2])
package_path <- as.character(args[3])
ichorCNA_outdir <- as.character(args[4])

# use current dir if no bamfile_path
if(is.null(bamfile_path) || is.na(bamfile_path)) {
  bamfile_path <- "."
}

# setup the output dir if not specified
if(is.null(ichorCNA_outdir) || is.na(ichorCNA_outdir)) {
	ichorCNA_outdir <- file.path(bamfile_path, "ichorCNA_outdir")
}

# create the output dir if not exist
if(!dir.exists(ichorCNA_outdir)){
	dir.create(ichorCNA_outdir, showWarnings = TRUE, recursive = TRUE)
	}

# setup the slurm partition if not specified
if(is.null(slurm_partition) || is.na(slurm_partition)) {
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

if(is.null(package_path) || is.na(package_path)) {
  package_path <- as.character("/home/nrlab/wang04/ulyses/explore_packaged")
  
}


# get the files

# file name pattern

# pattern1 is for the bams
pattern1 <- ".*0\\.1x\\.mrkdup\\.bam$"

# pattern2 is for the 1M bam files converted from finaledb files
pattern2 <- ".*hg19\\.frag\\.tsv\\.bgz\\.GRanges\\.rds\\.1M\\.rds\\.bam$"

bamfile1 <- list.files(path = bamfile_path,
	pattern = pattern1,
	recursive = TRUE,
	full.names = TRUE)

bamfile2 <- list.files(path = bamfile_path,
	pattern = pattern2,
	recursive = TRUE,
	full.names = TRUE)

bamfile <- c(bamfile1, bamfile2)
# prepare the rslurm input params for 0.2x bams

run_file_list <- tibble(`bamfile` = bamfile) |>
  dplyr::rename(`input_bam` = `bamfile` ) |>
  dplyr::mutate(wigfile = file.path(ichorCNA_outdir, paste0(basename(input_bam), ".wig"))) |>
  tibble::add_column(ichorCNA_outdir = ichorCNA_outdir)  |>
  dplyr::mutate(bam_basename = basename(input_bam)) |> 
  dplyr::mutate(param_file = paste0(ichorCNA_outdir, "/", bam_basename, ".ichor.params.txt")) |>
  dplyr::mutate(ichorCNA_script = case_when(
	stringr::str_detect(input_bam, pattern = "0\\.1x\\.mrkdup\\.bam$") ~ "/home/nrlab/resources/ichorCNA/rubicon/ichorcna_rubicon_haichao_self_control.sh",
    stringr::str_detect(input_bam, pattern = "1M\\.rds\\.bam$") ~ "/home/nrlab/resources/ichorCNA/rubicon/ichorcna_rubicon_haichao_self_control_finaleDB.sh"	
  ))


# filter by existence of rds files
run_file_list_filtered <- dplyr::filter(run_file_list, !file.exists(param_file))

# rslurm function
source(file.path(package_path,"functions/run_ichorCNA.R"))

# run rslurm_apply on the filtered list
#run_file_list_filtered <- run_file_list_filtered[nrow(run_file_list_filtered), ]

sjob <- rslurm::slurm_apply(
	f = run_ichorcna_self_ctrl, 
	params = run_file_list_filtered |> dplyr::select(-param_file),
	pkgs = c("tidyverse", "rslurm"),
	rscript_path = "/home/nrlab/tools/anaconda3/envs/R4_2/bin/Rscript",
	job_array_task_limit = 30,
	nodes = 3000, 
        cpus_per_node = 1,
        processes_per_node = 1,
	jobname = "ichorCNA3",
	slurm_options = list(time = '1-0',  mem = '16G', partition = slurm_partition),
	submit = TRUE)


message("Mission Completed, Goodbye!")
