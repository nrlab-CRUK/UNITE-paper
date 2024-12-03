
sbatch -J "bwa1 h sd100" -o mem_%j.out -e mem_%j.err --time=10:00:00 --mem=128G \
--wrap="source ~/.bashrc; conda activate /home/nrlab/tools/anaconda3/envs/projectx; bwa mem -t 40 -I 167,100 /scratcha/nrlab/resources/bwa_index/current/ucsc.hg19.fasta SLX-13900.D708_D504.HVK3YBBXX.r_1.fq.gz  SLX-13900.D708_D504.HVK3YBBXX.r_2.fq.gz | samtools sort -o SLX-13900.D708_D504.HVK3YBBXX.bam"




conda activate projectx

for f in SLX-13900.D708_D504.HVK3YBBXX.bam
do
sbatch -J "mrkdup ${f}" -o dup_%j.out -e dup_%j.err --time=1-0 --mem=128G --wrap="source ~/.bashrc; conda activate projectx; \
picard MarkDuplicates \
--INPUT ${f} \
--ASSUME_SORT_ORDER coordinate \
--REMOVE_DUPLICATES false \
--METRICS_FILE ${f}.mrkdup_metrics.txt \
--OUTPUT ${f}.mrkdup.bam \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY LENIENT"

done



WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in SLX-13900.D708_D504.HVK3YBBXX.bam.mrkdup.bam
do 
sbatch --time=2-0 --job-name='d_ulyses' --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 5000000"
done



WRAPPER_FILE="/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/wrapper_func2.r"

for i in PGDX8828P_WGS.sorted_processed.bam
do 
sbatch --time=2-0 --job-name='d_ulyses' --mem=64G --wrap="source ~/.bashrc; conda activate R4_1; /home/nrlab/tools/anaconda3/envs/R4_1/bin/Rscript ${WRAPPER_FILE}  $(pwd)/${i} 5000000"
done






####DELFI downsample

library(tidyverse)
library(rslurm)

setwd("/scratcha/nrlab/wang04/external_dataset_download/delfi/delfi_ega_data/bamfiles")
bamlist <- list.files(path = getwd(), full.names=T, pattern = "\\.*\\.bam$")
working_dir <- "/scratchb/nrlab/wang04/bwa_issue/delfi"
file_list <-  tibble::tibble(file_fullname = bamlist)

out_file <- file.path(working_dir, "all_bam_file_fullname.csv")

write_csv(file_list, file = out_file, col_names = FALSE)




master_file <- out_file

output_path <-  working_dir
output_csv_file <- file.path(output_path, "master_with_bam_metrics.csv")
short_lower <- 90
short_upper <- 150
long_lower <- 151
long_upper <- 210

# modify the master meta csv

master_before_mrkdup <- file_list %>%
    add_column(short_lower = short_lower, 
                short_upper = short_upper,
                long_lower = long_lower,
                long_upper = long_upper)



f.mrkdup_name <- function(x, path =  output_path ) {
    
    return(x)
}



pars <- dplyr::select(master, file_fullname, short_lower, short_upper, long_lower, long_upper)

slurm_options<- list(time = '20:00:00',  mem = '100G')

count_sjob <- slurm_apply(
    master_function, pars, 
    jobname = 'countBam',
    nodes = 2000, 
    cpus_per_node = 1,
    slurm_options = slurm_options,
    global_objects = c("count_all_bam", 
                        "count_good_bam", 
                        "count_dup_bam", 
                        "count_unmapped_bam", 
                        "count_short_bam", 
                        "count_long_bam", 
                        "count_chrM_bam"),
    submit = TRUE 
)





get_job_status(count_sjob)$queue
get_job_status(count_sjob)$log
# get_job_status(count_sjob)
# cancel_slurm(count_sjob)

#########################################################################
#get the result 
#########################################################################

count_result <- get_slurm_out(count_sjob, outtype = 'table', wait = TRUE)

#########################################################################
#analyze the result 
#########################################################################

join <- inner_join(count_result, master, by = "file_fullname")

join$n_nucleotides_all <- as.numeric(join$n_nucleotides_all)
join$n_reads_all <- as.numeric(join$n_reads_all)

join$n_nucleotides_good <- as.numeric(join$n_nucleotides_good)
join$n_reads_good <- as.numeric(join$n_reads_good)

join$n_nucleotides_short <- as.numeric(join$n_nucleotides_short)
join$n_reads_short <- as.numeric(join$n_reads_short)

join$n_nucleotides_long <- as.numeric(join$n_nucleotides_long)
join$n_reads_long <- as.numeric(join$n_reads_long)

join$n_nucleotides_dup <- as.numeric(join$n_nucleotides_dup)
join$n_reads_dup <- as.numeric(join$n_reads_dup)

join$n_nucleotides_unmapped <- as.numeric(join$n_nucleotides_unmapped)
join$n_reads_unmapped <- as.numeric(join$n_reads_unmapped)

join$chrM_count <- as.numeric(join$chrM_count)


join2 <- join %>%
    mutate(depth = n_nucleotides_good/3300000000,
        long_fraction = n_reads_long/n_reads_good,
        short_fraction = n_reads_short/n_reads_good,
        long_short_ratio = n_reads_long/n_reads_short,
	cohort = 'delfi',
	raw_data = 'y')

        #long_fraction_scaled = scale(join$long_fraction),
        #short_fraction_scaled = scale(join$short_fraction)) 

# save the joined files
write_csv(x = join2, output_csv_file )


# downsample


master <- join2 %>%
  dplyr::filter(depth >= 0.1)

target_depth_list <- as.list(seq(0.1, 1, 0.1))
outdir_root <- working_dir


create_dir <- function(x, outdir_root) {
  
  outdir2 <- file.path(outdir_root, paste0(x,'x'))
  dir.create(outdir2, showWarnings = TRUE, recursive = TRUE)
  return(0)  
}

lapply(target_depth_list, create_dir, outdir_root = outdir_root)



downsample_meta_form <- function(target_depth, outdir_root) {

  outdir2 <- file.path(outdir_root, paste0(target_depth,'x'))

  downsample_meta <- master  %>% 
    dplyr::filter(depth >= target_depth & raw_data == "y") %>%
    dplyr::select(file_fullname, depth, n_reads_good, cohort)  %>% 
    dplyr::mutate(n_reads_target = target_depth * 3000000000 / 100,
                downsamp_probability = n_reads_target/n_reads_good,
                target_depth = target_depth,
                #file_fullname_target = file.path(outdir2, paste0(basename(file_fullname), ".pos_based_downsamp_", target_depth, "x.bam")),
                file_fullname_target = file.path(outdir2, paste0(basename(file_fullname), ".downsamp_", target_depth, "x.bam")),
                mrkdup_output = paste0(gsub(".bam$", "", file_fullname_target), ".mrkdup.bam"))
                #pos_downsamp_mrkdup_output = paste0(gsub(".bam$", "", file_fullname_target), ".mrkdup.bam"))
}

downsample_meta_bind <- lapply(target_depth_list, downsample_meta_form, outdir_root = outdir_root)  %>% 
  bind_rows()

write_csv(x = downsample_meta_bind, file = file.path(outdir_root, paste0("downsample_meta_", "bind.csv")))



# rename the cols for rslurm parmas usage (DownsampleSamGATK)
downsample_meta_rename <- downsample_meta_bind %>%
  dplyr::rename(input = file_fullname,
    output = file_fullname_target,
    probability = downsamp_probability)  %>% 
  dplyr::select(input, output, probability, depth, mrkdup_output) %>%
  dplyr::mutate(mrkdup_input = output)






gatk <- function(gatk_path=gatk_path, method, args, maxheap="256G"){
  
  gatk_path <- "/scratcha/nrlab/wang04/urine/np_scripts/gatk-4.2.0.0/gatk"

  args <- cbind(paste("--java-options", paste0('"',"-Xmx", maxheap,'"'), method),
                args)
  retcode <- system2(gatk_path, args, stdout=FALSE)
  if(retcode != 0){
    stop( paste("GATK command [ java", paste(args, collapse=" "), "] failed."))
  }
  retcode
}


markDuplicatesGATK <- function(input, 
  output, 
  create_index = TRUE, 
  metrics_file = paste0(output, ".metrics"),
  quiet = FALSE,
  stringency = "LENIENT",
  verbosity = "ERROR",
  overwrite_output = FALSE) {
                                
 if(!file.exists(input)) stop(paste("markDuplicatesGATK cannot find", input))

 if (!overwrite_output) {
    if(file.exists(output)) stop("markDuplicatesGATK failed to create", output, "file already exists.")
 } else {
    unlink(output)
 }

  
  args <- paste("--INPUT", input,
                "--OUTPUT", output,
                "--CREATE_INDEX", create_index,
                "--METRICS_FILE", metrics_file,
                "--VERBOSITY", verbosity,
                "--VALIDATION_STRINGENCY", stringency,
                "--QUIET", quiet)

  gatk(method = "MarkDuplicates", args = args)

}





downsampleSamGATK <- function(input, 
  output, 
  probability, 
  depth,
  strategy = "HighAccuracy", 
  accuracy = 0.00001, 
  create_index = TRUE, 
  metrics_file = paste0(output, ".metrics"),
  quiet = FALSE,
  stringency = "LENIENT",
  verbosity = "ERROR",
  overwrite_output = FALSE) {
                                
 if(!file.exists(input)) stop(paste("downsampleSamGATK cannot find", input))

 if (!overwrite_output) {
    if(file.exists(output)) stop("downsampleSamGATK failed to create", output, "file already exists.")
 } else {
    unlink(output)
 }

 if (as.numeric(depth) > 10) strategy <- "Chained"
  
  args <- paste("--INPUT", input,
                "--OUTPUT", output,
                "--STRATEGY", strategy,
                "--PROBABILITY", probability,
                "--ACCURACY", accuracy,
                "--CREATE_INDEX", create_index,
                "--METRICS_FILE", metrics_file,
                "--VALIDATION_STRINGENCY", stringency,
                "--VERBOSITY", verbosity,
                "--QUIET", quiet)

  gatk(method = "DownsampleSam", args = args)

}




downsampleSamGATK_markDuplicatesGATK <- function(input, 
  output, 
  probability, 
  depth,
  mrkdup_input, 
  mrkdup_output) {
  
 downsampleSamGATK(input = input, output = output, probability = probability, depth = depth)
 markDuplicatesGATK(input = mrkdup_input, output = mrkdup_output)
 unlink(mrkdup_input)
 unlink(paste0(gsub(".bam", "", mrkdup_input), ".bai"))

}



sjob2 <- slurm_apply(f = downsampleSamGATK_markDuplicatesGATK, 
    params = downsample_meta_rename, 
    jobname = 'downsampGATK',
    nodes = 2000, 
    cpus_per_node = 1,
    slurm_options = list(time = '20:00:00',  mem = '100G'),
    global_objects = c("gatk", "markDuplicatesGATK", "downsampleSamGATK"),
    submit = TRUE)


