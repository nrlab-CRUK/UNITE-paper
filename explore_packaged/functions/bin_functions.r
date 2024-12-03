

set_bin_index <- function(combine_n_bins, object, armlevels){
  
  # for p and q arms, count from telomere to centromere

  r1 <- object %>%
    plyranges::group_by(arm) %>%
    plyranges::mutate(arm_bin_index = 1:length(arm)) %>%
    plyranges::mutate(arm_bin_index = ifelse(stringr::str_detect(arm, "p"), 
                                              ceiling(arm_bin_index/!!combine_n_bins),
                                              ceiling(rev(arm_bin_index/!!combine_n_bins)))) %>%
    plyranges::ungroup()
  
  # set levels of arm 
  
  r1$arm <- factor(r1$arm, levels = armlevels)
  
  r1 <- r1 %>%
    plyranges::mutate(id = paste0(arm, "_", arm_bin_index))
  

  
  # remove redundant bins
  
  discard_bin <-  r1 %>%
    plyranges::group_by(arm, arm_bin_index) %>%
    plyranges::summarise(count = n()) %>%
    tibble::as_tibble() %>%
    dplyr::filter(count < !!combine_n_bins) %>%
    dplyr::mutate(id  = paste0(arm, "_", arm_bin_index))
  
  # remove incomplete bins from r1
  
  r2 <- r1 %>%
    plyranges::filter(!id %in% discard_bin$id)
  
  # correct the order/direction of 
  
  r3 <- r2 %>%
    plyranges::group_by(arm) %>%
    plyranges::mutate(arm_bin_index  = ifelse(stringr::str_detect(arm, "p"), arm_bin_index, rev(arm_bin_index))) %>%
    plyranges::mutate(id = paste0(arm, "_", arm_bin_index)) %>%
    plyranges::ungroup()
  
  # set the levels of arm_bin_index and id
  
  r3$arm_bin_index <- factor(r3$arm_bin_index, levels = sort(unique(r3$arm_bin_index)))
  
  id_levels <- tibble::as_tibble(r3) %>%
    dplyr::select(arm, arm_bin_index, id) %>%
    dplyr::arrange(arm, arm_bin_index) %>%
    dplyr::select(id) %>%
    dplyr::pull(id) %>%
    unique()
  
  r3$id <- factor(r3$id, levels = id_levels)
  
  

  return(r3)
}

