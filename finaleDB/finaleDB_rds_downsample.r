

library(rslurm)
library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)

###############################################################################
# params
###############################################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least rdsfile_path must be supplied (input file).\n", call.=FALSE)
}

# rdsfile
rdsfile_path <-  as.character(args[1])
n_frag_selected <- as.numeric(args[2])
outdir <-  as.character(args[3])
slurm_partition <- as.character(args[4])

# get the files

rdsfile <- list.files(path = rdsfile_path, 
	pattern = "*.frag.tsv.bgz.GRanges.rds$", 
	full.names = TRUE)


run_file_list <- tibble(rdsfile = rdsfile, n_frag_selected = n_frag_selected ) %>%
	dplyr::mutate(outfile = file.path(outdir, paste(basename(rdsfile), ".", n_frag_selected, "M", ".rds", sep = "")))

# filter by existence of rds files

run_file_list_filtered <- dplyr::filter(run_file_list, 
	!file.exists(outfile))

# helper function

rds_downsample_helper <- function(rdsfile, n_frag_selected, outfile){

	n_frag_selected_true <- as.numeric(n_frag_selected) * 1000000
	frag <- readRDS(rdsfile)

	frag

	if(length(frag) >= n_frag_selected_true){
		ans <- sample(frag, n_frag_selected_true, replace = FALSE)
		message("Saving to ", outfile)
		saveRDS(ans, outfile)
	} else {
		stop("n_frag_selected is larger than the number of fragments in the rds file.")
	}

}


# run rslurm_apply on the filtered list

sjob <- rslurm::slurm_apply(
	f = rds_downsample_helper, 
	params = run_file_list_filtered,
	pkgs = c("tidyverse", "rslurm", "GenomicRanges", "GenomicAlignments"),
	rscript_path = "/home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript",
	job_array_task_limit = 300,
	nodes = 2000, 
        cpus_per_node = 1,
        processes_per_node = 1,
	jobname = "rds_downsamp_5d",
	slurm_options = list(time = '5-0',  mem = '32G' , partition = slurm_partition),
	submit = TRUE)


message("Mission Completed, Goodbye!")