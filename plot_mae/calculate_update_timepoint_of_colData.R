
library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)
library(argparse)

# make argument parser
parser <- argparse::ArgumentParser()
# add input_mae_file argument, set default to MAE4.rds
parser$add_argument("--input_mae_file", default = "/scratchb/nrlab/wang04/ulyses_feat_viz/MAE4.rds")
# add output_mae_file argument, set default to MAE5.rds
parser$add_argument("--output_mae_file", default = "/scratchb/nrlab/wang04/ulyses_feat_viz/MAE5.rds")
# add updated_colData argument, set default to /home/nrlab/wang04/ulyses/meta_data/colData.csv
parser$add_argument("--updated_colData", default = "/home/nrlab/wang04/ulyses/meta_data/colData.csv")
# add mae_tier2_path argument, set default to /mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/
parser$add_argument("--mae_tier2_path", default = "/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/")


# parse arguments
args <- parser$parse_args()
input_mae_file <- args$input_mae_file
output_mae_file <- args$output_mae_file
updated_colData <- args$updated_colData
mae_tier2_path <- args$mae_tier2_path


# data

message("Loading MAE RDS data")
mae <- readRDS(input_mae_file)
source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

# get the colData from input_mae_file
mae_coldata <- colData(mae) %>%
  as.data.frame() %>%
  # only keep these columns: tmad
  select(c("tmad")) 
mae_rownames <- mae_coldata %>% 
  rownames()


# get the updated colData
updated_coldata <- read.csv(updated_colData)
# remove bam_id col is.na
updated_coldata <- updated_coldata %>%
  filter(!is.na(bam_id)) %>%
  # convert bam_id to rownames
  column_to_rownames(var = "bam_id")



# merge the updated colData with the mae_coldata
coldata_new <- merge(mae_coldata, updated_coldata, by = "row.names", sort = FALSE)

coldata_new$age <- as.numeric(coldata_new$age)
coldata_new$clinical_tf <- as.numeric(coldata_new$clinical_tf)

rownames(coldata_new) <- coldata_new$Row.names
coldata_new$Row.names <- NULL

# rebuild MAE3 with new colData
colData(mae) <- as(coldata_new, "DataFrame")

# save MAE3 with new colData as MAE4.rds
mae_updated <- mae
saveRDS(mae_updated, output_mae_file)

# backup the new mae file 
# copy output_mae_file to mae_tier2_path
file.copy(output_mae_file, mae_tier2_path)