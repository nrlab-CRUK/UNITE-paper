## download, trim and align using TAP pipeline

```bash

 /home/bioinformatics/software/pipelines/kickstart/current/bin/kickstart   -f "Library Type" -f "Index Type" \
 -l=SLX-20706 \
 -l=SLX-22151 \
 -l=SLX-21409 \
 -l=SLX-22439 \
 -l=SLX-22440 \
 -l=SLX-19894 \
 -l=SLX-21959 \
 -l=SLX-21619

```

## calculate bam stats

```
cd /mnt/scratchb/nrlab/wang04/dbs/processed

for i in *.bam

do

sbatch --time=0-10 --mem=128G  -J "count" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/summariseBam.R  ${i}  .  FALSE"

done

```

## downsample to 0.1x (1Million fragments)

```bash

cd /mnt/scratchb/nrlab/wang04/dbs/processed
downsampled_dir=/mnt/scratchb/nrlab/wang04/dbs/downsampled_bam

for i in *.csv

do

sbatch --time=0-20 --mem=32G  -J "downsamp" -p "general" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/downsample_gatk.R ${i}  0.1  ${downsampled_dir}  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"

done



```


## run ulyses: generate ulyses_image.rds files

```bash
downsampled_dir=/scratchc/nrlab/wang04/ulyses/data/internal/dbs/downsampled_bam

sbatch --time=0-00:30:00  --job-name='u-w-igniter' --partition="epyc" --mem=1G --wrap="source ~/.bashrc; conda activate R4_1; Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ulyses.R  ${downsampled_dir}"
```

## ulyses plot

```bash
cd /scratchc/nrlab/wang04/ulyses/data/internal/dbs/downsampled_bam

conda activate R4_2

for i in *.ulyses_image.rds
do 
Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R   $i
done

```

## run ichorCNA


```bash

# DBS sample ichor

conda activate R4_2

Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ichorCNA.R /scratchc/nrlab/wang04/ulyses/data/internal/dbs/downsampled_bam  epyc  /home/nrlab/wang04/ulyses/explore_packaged /scratchc/nrlab/wang04/ulyses/data/internal/dbs/ichorCNA_outdir2

```

## build MAE


### meta data
```R
library(tidyverse)
#ulyses_package_dir <- "/Users/wang04/Documents/libs/github/ulyses"
ulyses_package_dir <- "/home/nrlab/wang04/ulyses"
#out_dir <- "/Users/wang04/Documents/libs/github/ulyses/meta_data"
out_dir <- file.path(ulyses_package_dir, "meta_data", "DBS")

colData_cols <- c("cohort",
                  "patient_id",
                  "specimen_id",
                  "timepoint",
                  "bam_id",
                  "author",
                  "sample_type",
                  "ichorcna_tf",
                  "ichorcna_ploidy",
                  "ichorcna_gender",
                  "stage",
                  "library_kit",
                  "extraction_kit",
                  "seq_platform",
                  "age",
                  "gender",
                  "seq_depth",
                  "clinical_tf",
                  "data_source",
                  "mito_frac"
                  )

add_missing_cols <- function(cohort_colData, colData_cols) {
  
  tmp <- setdiff(colData_cols, colnames(cohort_colData)) 
  tmp_named <- setNames(tmp, tmp)
  cols_to_add <- dplyr::bind_rows(tmp_named)[0, ]
  cols_to_add[seq_len(nrow(cohort_colData)), ] <- NA
  print("Missing cols:")
  print(tmp)
  
  ans <- dplyr::mutate(cohort_colData, cols_to_add) %>%
    dplyr::select(all_of(colData_cols))
  return(ans)
}

######################################################
# read meta data
######################################################

dbs_meta_file <- file.path(ulyses_package_dir, "meta_data", "DBS", "final", "dbs_meta_mannual.csv")

ulyses_colData <- read_csv(dbs_meta_file)


######################################################
# read coverage information
######################################################

coverage_path <- "/scratchc/nrlab/wang04/ulyses/data/internal/dbs/"

depth_fl <- list.files(path = coverage_path , pattern = ".bam.summary.csv", recursive=TRUE, full.names = TRUE)

depth_dt <- purrr::map_df(depth_fl, read_csv, show_col_types = FALSE)

seq_depth_dt <- depth_dt %>%
  mutate(bam_id = str_extract(file, "(^[^.]+\\.[^.]+)")) %>%
  select(all_of(c("bam_id", "n_read", "n_read_mapped", "n_read_duplicate", "n_read_mapped_chrM"))) %>%
  mutate(n_read_mapped_unique = n_read_mapped - n_read_duplicate) %>% 
  mutate(n_frag_unique = floor(n_read_mapped_unique / 2)) %>%
  mutate(seq_depth = n_frag_unique / 1E7) %>%
  mutate(mito_frac = n_read_mapped_chrM / n_read_mapped) %>%
  select(all_of(c("bam_id", "seq_depth", "mito_frac")))


ulyses_colData <- left_join(ulyses_colData, 
  seq_depth_dt, by = "bam_id")


######################################################
# add ichorCNA_tf information
######################################################

ichorcna_path <- "/scratchc/nrlab/wang04/ulyses/data/internal/dbs/ichorCNA_outdir"

read_param <- function(x) {
    if(file.exists(x)) {
    #file.copy(x, "./tmp")
    a <- read.delim(x, header = TRUE, nrows = 1) %>%
        as_tibble()

    gender <- readLines(x, n = 6, warn = FALSE)[[6]] %>%
    stringr::str_remove(pattern = "^Gender:\t")

    ans <- a %>%
    add_column(`gender` = gender)
    return(ans)

    } else {
        message(paste(x, "does not exist!"))
    }

}



ichorcna_fl <- list.files(path = ichorcna_path, pattern = ".ichor.params.txt", recursive = TRUE, full.names = TRUE)

ichorcna_tf_dt <- purrr::map_df(ichorcna_fl, read_param )



ichorcna_tf_dt_tidy <- ichorcna_tf_dt %>%
  dplyr::mutate(bam_basename = str_remove(Sample, ".ichor$")) %>% 
  dplyr::rename(`ichorcna_tf` = `Tumor.Fraction`) %>% 
  dplyr::rename(`ichorcna_ploidy` = `Ploidy`) %>%
  dplyr::rename(`ichorcna_gender` = `gender`) %>% 
  dplyr::mutate(bam_id = str_extract(bam_basename, "(^(SLX|DL)\\-\\d+\\.[^.]+)|(^EE\\d+)")) %>%
  dplyr::select(all_of(c("bam_id", "ichorcna_tf", "ichorcna_ploidy", "ichorcna_gender"))) 

ulyses_colData <- left_join(ulyses_colData, 
  ichorcna_tf_dt_tidy, by = "bam_id")


######################################################
# add missing cols 
######################################################

ulyses_colData <- add_missing_cols(cohort_colData = ulyses_colData,
                 colData_cols = colData_cols)


######################################################
# save the results
######################################################
saveRDS(ulyses_colData, file = file.path(out_dir,"dbs_colData.rds"))
write_csv(ulyses_colData, file = file.path(out_dir, "dbs_colData.csv"))


```


## build MAE

```bash

dbs_wd=/scratchc/nrlab/wang04/ulyses/data/internal/dbs
cd $dbs_wd
cp /home/nrlab/wang04/ulyses/multiassayexperiment/ulyses_mae.R $dbs_wd


screen -RD dbs_mae

sintr --time=3-0 --mem=200G -J 'all_mae' --partition="epyc"

conda activate R4_2

Rscript ulyses_mae.R

```

## prepare the test tensor

```bash
# get the data

wd=/Users/wang04/cnn_all_dbs_tp1
mkdir -p ${wd}
cd ${wd}
cp /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R ${wd}

# all timepoints
mae_file=/scratchc/nrlab/wang04/ulyses/data/internal/dbs/dbs_MAE.rds

sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript  ./scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1 --sample_type fpDBS  --outdir ${wd}  --mae_file ${mae_file} "

# tp1

mae_file=/scratchc/nrlab/wang04/ulyses/data/internal/dbs/dbs_MAE.rds

sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript  ./scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1 --sample_type fpDBS  --outdir ${wd}  --mae_file ${mae_file} --timepoints 1  --authors 'A.A. et al' "

==============================

```


## apply the FM_SC_tp1 model 

```bash

# train
train_wd=/Users/wang04/cnn_all_FM_SC_tp1
test_wd=/Users/wang04/cnn_all_dbs_tp1
cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}

cd $train_wd


# simpler is better in PP tp1 test
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 35 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type simpler --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"





# test


cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.simpler.epochs_35.batch_size_32.verbose_1.validation_split_0.4.learning_rate_0.001.h5"


# test using PP_tp1

test_wd=/Users/wang04/cnn_all_PP_tp1
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.simpler.epochs_35.batch_size_32.verbose_1.validation_split_0.4.learning_rate_0.001.h5"


```





