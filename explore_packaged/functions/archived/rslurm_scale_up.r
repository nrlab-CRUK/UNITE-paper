library(tidyverse)

setwd("/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/bams")

wrapper_file <- "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func.r"
meta_table_file <- "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/golden_master.csv"
source(wrapper_file)

# prepare rslurm file list



meta <- read_csv(meta_table_file, show_col_types = FALSE) %>% 
        dplyr::filter(target_depth == 0.5) %>% 
        dplyr::filter(specimen == "urine") %>% 
        dplyr::select(mrkdup_output) %>%
        dplyr::rename(bamfile = mrkdup_output)


global_objects <- ls() %>% as.vector()

# run
sjob_scale_up2 <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta,
                  jobname = 'upscale2',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '30-0',  mem = '8G'),
                  job_array_task_limit = 40,
                  global_objects = global_objects,
                  submit = TRUE)

# check the running results
 rslurm::get_job_status(sjob_scale_up2)

# run
sjob_scale_up3 <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta,
                  jobname = 'upscale3',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '30-0',  mem = '8G'),
                  job_array_task_limit = 40,
                  global_objects = global_objects,
                  submit = TRUE)

# check the running results
 rslurm::get_job_status(sjob_scale_up3)


 # run
sjob_scale_up4 <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta,
                  jobname = 'upscale4',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '30-0',  mem = '8G'),
                  job_array_task_limit = 40,
                  global_objects = global_objects,
                  submit = TRUE)

sjob_scale_up5 <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta,
                  jobname = 'upscale5',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '30-0',  mem = '8G'),
                  job_array_task_limit = 40,
                  global_objects = global_objects,
                  submit = TRUE)


sjob_scale_up6 <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta,
                  jobname = 'upscale6',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '1-0',  mem = '4G'),
                  job_array_task_limit = 100,
                  global_objects = global_objects,
                  submit = TRUE)
sjob_scale_up7 <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta,
                  jobname = 'upscale7',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '1-0',  mem = '4G'),
                  job_array_task_limit = 10,
                  global_objects = global_objects,
                  submit = TRUE)


# check the running results
 rslurm::get_job_status(sjob_scale_up4)


 # run a plasma sample
meta_plasma <- tibble("bamfile" = "/scratcha/nrlab/wang04/plasma/SLX-11379.D704_D508.trimgalore.bwamem2_mrkdup_filtered.bam")

sjob_scale_up_plasma <- rslurm::slurm_apply(
                  f = bam_to_image,
                  meta_plasma,
                  jobname = 'up_plasma',
                  nodes = 2000, 
                  cpus_per_node = 1,
                  slurm_options = list(time = '30-0',  mem = '8G'),
                  job_array_task_limit = 40,
                  global_objects = global_objects,
                  submit = TRUE)
