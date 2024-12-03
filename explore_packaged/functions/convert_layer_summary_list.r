

library(tidyverse)
library(reticulate)
file <- "/scratchc/nrlab/wang04/ulyses_second_batch/try_ulyses_gpu/SLX-18445.SXTLI192.HiSeq4000.bam.0.1x.mrkdup.bam.ulyses_bind_olap_chunks.rds.pkl"

# set the default reticulate python env

use_condaenv("R4_2")
pd <- import("pandas")
summary_list <- pd$read_pickle(file)


export_summary_table_long <- function(x, 
                                      layers){
  
  
  tibble_wide <-  dplyr::bind_rows(x, .id = "tile_size") %>%
    dplyr::mutate(unique_group_id = row_number()) %>%
    dplyr::mutate(chr = stringr::str_extract(id, pattern = "\\b\\d{1,2}") %>% as.integer()) %>%
    dplyr::mutate(tile_size_label = stringr::str_extract(tile_size, pattern = "\\d+\\b") %>% as.integer()) %>%
    tidyr::separate(col = id, into = c("arm", "within_arm_index"), sep = "_", remove = FALSE, convert = TRUE) %>%
    dplyr::relocate(chr, .after = arm) 
  
  tibble_wide$arm <- factor(tibble_wide$arm, levels = tibble_wide$arm %>% unique()) 
  tile_size_levels <- paste0("n_100kb_bins_", tibble_wide %>% dplyr::pull(tile_size_label) %>% sort() %>% unique() )  
  tibble_wide$tile_size <-  factor(tibble_wide$tile_size, levels = tile_size_levels)
  
  # arrange columns in specific order
  tibble_wide <- tibble_wide %>%
    dplyr::arrange(tile_size, id, within_arm_index, frag_len) %>%
    # careful
    dplyr::group_by(tile_size, frag_len) %>% 
    dplyr::mutate(id_rename = row_number() ) %>% 
    dplyr::relocate(id_rename, .after = id) %>%
    dplyr::relocate(tile_size_label, .after = tile_size) %>%
    dplyr::ungroup()
  
  tibble_long <- tibble_wide %>%
    tidyr::pivot_longer( 
      cols = all_of(layers), 
      names_to = "channel", 
      values_to = "pixel") %>%
    dplyr::mutate(unique_group_id = row_number())
  
  tibble_long$channel <- factor(tibble_long$channel , levels = layers)
  
  tibble_long <- tibble_long %>%
    # the order of cols in the arrange correspond to the 
    # dimension of ndarray in python analysis
    dplyr::arrange(tile_size, channel, id, frag_len)

  return(tibble_long)
  
}



tibble_long <- export_summary_table_long(x = summary_list, 
                                           layers = layers)
  
# calculate ulyses object
tibble_long_grouped <- tibble_long %>% 
    # the order is important
    dplyr::arrange(tile_size, channel,frag_len, id) %>%
    dplyr::group_by(tile_size)
  
  