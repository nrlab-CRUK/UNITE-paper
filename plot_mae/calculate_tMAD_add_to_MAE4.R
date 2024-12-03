
library(tidyverse)
library(patchwork)
library(MultiAssayExperiment)
library(nord)
library(plotly)
library(htmlwidgets)
library(ggpubr)

# function
tmad_seg <- function(x){
  CNA.object <- DNAcopy::CNA(genomdat = x$value, 
                             chrom = x$chr,
                             maploc = x$pos,
                             data.type = "logratio",
                             sampleid = unique(x$primary),
                             presorted = TRUE)
  
  smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object)
  
  segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, 
                                                  undo.splits = "sdundo",
                                                  undo.SD = 3,
                                                  verbose=1)
  
  seg_num.mark <- segment.smoothed.CNA.object$output$num.mark
  seg_seg.mean<- segment.smoothed.CNA.object$output$seg.mean
  
  seg_value <- rep(seg_seg.mean, seg_num.mark)
  
  ans <-  dplyr::mutate(x, seg = seg_value)
  return(ans)
}


tmad_function <- function(x){ 
  abs(mad(x = x, center=0, na.rm = TRUE))
}



# data



message("Loading MAE RDS data")
mae <- readRDS("/scratchb/nrlab/wang04/ulyses_feat_viz/MAE3.rds")
source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")
noncancer_disease <- noncancer_disease_hyper
outliers <- outliers_hyper



assay_to_keep <- c("cnv")

message("Filtering samples")
all <- mae[, , assay_to_keep]



all_long <- longFormat(all, 
            colDataCols = c("cohort", 
                            "author", 
                            "timepoint",
                            "clinical_tf",
                            "ichorcna_tf"
                            )) %>%
  tibble::as_tibble() 

all_long_filter <- all_long 



# add seg value to tibble
all_cnv_seg <- all_long_filter %>%
  dplyr::filter(assay == "cnv") %>%
  dplyr::filter(rowname != "6p_6" ) %>%
  dplyr::mutate(chr = str_extract(rowname, "\\d+")) %>%
  group_by(primary, chr) %>%
  dplyr::mutate(pos = row_number()) %>%
  ungroup() %>%
  dplyr::group_by(primary)

all_cnv_seg$rowname <- factor(all_cnv_seg$rowname, levels = all_cnv_seg$rowname %>% unique())

all_cnv_seg_list <- group_split(all_cnv_seg)
names(all_cnv_seg_list) <- group_keys(all_cnv_seg) %>% pull()

all_cnv_seg <- lapply(all_cnv_seg_list, tmad_seg)

ans <- all_cnv_seg %>%
  bind_rows() %>%
  select(-c(chr, pos))


# calculate tMAD
ans_tmad <- group_by(ans, primary) %>%
  mutate(tmad = tmad_function(seg)) 

triple <- ans_tmad %>% 
  group_by(primary) %>%
  summarise(tmad = unique(tmad), 
            ichorcna_tf = unique(ichorcna_tf), 
            clinical_tf = unique(clinical_tf))


# remove rows with NA usine filter across all columns
triple_filtered <- triple %>% 
  dplyr::filter_all(all_vars(!is.na(.)))


triple_filtered$clinical_tf <- as.numeric(triple_filtered$clinical_tf)
# pivot_longer 
triple_filtered_long <- triple_filtered %>%
  pivot_longer(cols = c("tmad", "ichorcna_tf"), 
               names_to = "metrics", 
               values_to = "predicted") %>%
  mutate(variable = factor(metrics, levels = c("tmad", "ichorcna_tf")))


# get the colData from MAE3 and add seg and tMAD to the colData

to_merge <- all_tmad_ichor <- triple %>%
  dplyr::select(primary, tmad) %>%
  column_to_rownames(var = "primary")

to_merge_rownames <- rownames(to_merge)


mae_rownames <-colData(mae) %>%
  as.data.frame() %>% 
  rownames()

mae_coldata <- colData(mae) %>%
  as.data.frame()

coldata_new <- merge( mae_coldata, to_merge, by = "row.names", sort = FALSE)

coldata_new$age <- as.numeric(coldata_new$age)
coldata_new$clinical_tf <- as.numeric(coldata_new$clinical_tf)

rownames(coldata_new) <- coldata_new$Row.names
coldata_new$Row.names <- NULL

# rebuild MAE3 with new colData
colData(mae) <- as(coldata_new, "DataFrame")

# save MAE3 with new colData as MAE4.rds
mae4 <- mae

saveRDS(mae4, "/scratchb/nrlab/wang04/ulyses_feat_viz/MAE4.rds")

# copy MAE4.rds to "/mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets"
system("cp /scratchb/nrlab/wang04/ulyses_feat_viz/MAE4.rds /mnt/nas-data/nrlab/group_folders/wang04/ulyses_datasets/")



# add timepoint information to the K.S and P.J. et al dataset
