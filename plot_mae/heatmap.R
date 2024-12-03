
library(ComplexHeatmap)

library(tidyverse)

# handling big dataframe --------------------------------------------------
heatmap_annotation_map <- all_long_filter %>% 
  filter(assay == "sd") %>%
  group_by(primary, 
         ichorcna_tf,
         ichorcna_tf_strat,
         cohort, 
         author, 
         stage, 
         library_kit, 
         extraction_kit, 
         seq_platform,
         age,
         gender,
         data_source,
         bicohort
         ) %>%
  summarise(n = n()) %>%
  select(-n)


# make age col numeric
heatmap_annotation_map$age <- as.numeric(heatmap_annotation_map$age)

# homogenize gender 
heatmap_annotation_map <- heatmap_annotation_map %>%
  mutate(gender = case_when(
    gender == "female" ~ "F",
    gender == "male" ~ "M",
    TRUE~gender
  ))
# homogenize seq_platform
heatmap_annotation_map <- heatmap_annotation_map %>%
  mutate(seq_platform = case_when(
    str_detect(seq_platform, "Novaseq_\\d$") ~ "NovaSeq",
    str_detect(seq_platform, "HighSeq_\\d$") ~ "HiSeq",
    str_detect(seq_platform, "HiSeq \\d+$") ~ "HiSeq",
    str_detect(seq_platform, 'HiSeq 2000/2500') ~ "HiSeq",
    str_detect(seq_platform, 'Novaseq_2 + HighSeq_1') ~ "NovaSeq/HiSeq",
    TRUE~seq_platform
  ))

# make <chr> column <fct>
heatmap_annotation_map <- heatmap_annotation_map %>% 
  ungroup() %>%
  mutate_if(is.character, as.factor) %>%
  group_by(ichorcna_tf_strat)%>%
  arrange(ichorcna_tf_strat, ichorcna_tf)

# set levels 

heatmap_annotation_map$bicohort <- factor(heatmap_annotation_map$bicohort, levels = c("Healthy", "Cancer"))


# rename columns

heatmap_annotation_map <- heatmap_annotation_map %>%
  rename('Cohort' = 'cohort',
         'TF Strat' = 'ichorcna_tf_strat',
         'Author' = 'author',
         'Stage' = 'stage',
         'Library Kit' = 'library_kit',
         'Extraction Kit' = 'extraction_kit',
         'Sequencer' = 'seq_platform',
         'Age' = 'age',
         'Gender' = 'gender',
         'Data Source' = 'data_source',
         'Binary Cohort' = 'bicohort')






get_tf <- function(x){
  
  tf_long_tibble2 <- select(tf_long_tibble, primary, ichorcna_tf) %>%
    group_by(primary) %>%
    summarise(ichorcna_tf = unique(ichorcna_tf))
  
  if(x %in% tf_long_tibble2$primary){
    
  ans <- tf_long_tibble2[tf_long_tibble2$primary == x,]$ichorcna_tf
    
  } else {
    ans <- NA_real_
  }
  return(ans)
}

plot_heatmap <- function(input, 
                         index, 
                         annotation_name_side = "right",
                         show_heatmap_legend = TRUE, 
                         legend_title = "z-score",
                         legend_at = seq(-5, 5, 2),
                         show_chr_name = TRUE,
                         chr_anno_text_size = 7,
                         show_detailed_sample_annotation = TRUE,
                         show_annotation_name = TRUE,
                         show_point_anno_axis = TRUE){
  
  input_mat <- input[[index]]
  
  
 ###############################
 # make left annotation
 ###############################

  panel_fun = function(index, nm) {
    #grid.rect()
    grid.segments(x0 = unit(1, "npc"), 
                  y0 = unit(0, "npc"),
                  x1 = unit(1, "npc"), 
                  y1 = unit(1, "npc"),
                  gp = gpar(col = "darkgrey")
                  )
    
    grid.segments(x0 = unit(0.8, "npc"), 
                  y0 = unit(0, "npc"),
                  x1 = unit(1, "npc"), 
                  y1 = unit(0, "npc"),
                  gp = gpar(col = "darkgrey")
                  )
    
    grid.segments(x0 = unit(0.8, "npc"), 
                  y0 = unit(1, "npc"),
                  x1 = unit(1, "npc"), 
                  y1 = unit(1, "npc"),
                  gp = gpar(col = "darkgrey")
                  )
    
    grid.text(nm, 
              x = 0.28, 
              y = 0.5, 
              gp=gpar(col="black", 
                                    fontsize= chr_anno_text_size, 
                                    fontface = "bold"))
    
  }
 
  # make align_to index
  
  
  rn <- input_mat |> rownames() |> stringr::str_extract("\\d+")
  rn_tibble <- tibble(chr = rn) %>%
    mutate(index = row_number()) %>%
    mutate(chr = paste0("Chr ", chr)) %>%
    group_by(chr)
  
  
  ans <- rn_tibble %>% 
    group_by(chr) %>% 
    summarise(min = min(index), max = max(index)) %>%
    mutate(range = paste0(min, ":", max))
  
  
  tmp <- rn_tibble %>% 
    group_split()
  
  names(tmp) <- group_keys(rn_tibble) %>% pull()
  
  align_to <- lapply(tmp, pull, 2)
  
 chr_name <- rowAnnotation(CHR = anno_block(
   align_to = align_to,
   panel_fun = panel_fun,
   width = unit(0.6, "cm")
 ),
 border = FALSE)

 if(show_chr_name == FALSE) {
   chr_name <- NULL
 }
 
 ###############################
 # make top annotation
 ###############################
 #ichorcna_tf <- sapply(colnames(input_mat), get_tf)
 
 
 # order ichorcna_tf
 #ichorcna_tf_order <- ichorcna_tf[order(ichorcna_tf)]
 # reorder input_mat based on ichorcna_tf
 #input_mat_order <- input_mat[, order(ichorcna_tf)]
 
 #col_boolean <- !is.na(ichorcna_tf_order)
 #ichorcna_tf_order <- ichorcna_tf_order[col_boolean]
 #input_mat_order <- input_mat_order[, col_boolean]
 
 
 boo <- colnames(input_mat) %in% unique(heatmap_annotation_map$primary)
 
 input_mat2 <- input_mat[, boo]
 
 selected <- heatmap_annotation_map %>%
   ungroup() %>%
   filter(primary %in% colnames(input_mat2))  
 
 col_order <- pull(selected, var = primary) 
 col_order <- factor(col_order, levels = col_order) %>%
   as.vector()
 
 input_mat_order <- input_mat2[,col_order]
 ichorcna_tf_order <- pull(selected, var = ichorcna_tf)
 
 
 
#  row_title <- paste(dim(input_mat_order)[[1]], " ", "5-Mb bins", sep = "")
 row_title <- paste( "5 Mb Genomic Bins", "(n=", dim(input_mat_order)[[1]], ")", sep = "")
 

 col_title <- paste(names(input)[[index]], "\n (", "n =", dim(input_mat_order)[[2]], ")")
 
 tf_simple_anno_axis_param = list(
   side = "left",
   at = c(0, 0.8),
   facing = "inside",
   gp = gpar(fontsize = 5)
 )
 
 
 tf_anno <-  anno_points(ichorcna_tf_order,
                           border = FALSE,
                           gp = gpar(col = "orange3"),
                           size = unit(1, "mm"),
                           ylim = c(0, 0.8),
                          extend = 0.01,
                         axis = show_point_anno_axis,
                           height = unit(0.5, "cm"),
                           axis_param = tf_simple_anno_axis_param
                           )
 
 
 
 # make mad anno
 
 mad_simple_anno_axis_param = list(
   side = "left",
   at = c(0, 0.5, 1.0),
   facing = "inside",
   gp = gpar(fontsize = 5)
 )
tMAD_func <- function(x){ 
  abs(mad(x = x, center=0, na.rm = TRUE))
}

 mad_anno <-  anno_points(apply(input_mat_order, FUN=tMAD_func, MARGIN = 2),
                           border = FALSE,
                           gp = gpar(col = "grey"),
                           size = unit(1, "mm"),
                          extend = 0.01,
                         axis = show_point_anno_axis,
                           ylim = c(0, 1.0),
                           height = unit(0.5, "cm"),
                           axis_param = mad_simple_anno_axis_param
                           )
 
 # order the df anno 
# heatmap_annotation_map_selected <- heatmap_annotation_map %>%
#   dplyr::filter(primary %in% colnames(input_mat_order))
 
# heatmap_annotation_map_selected$primary <- factor(heatmap_annotation_map_selected$primary, 
#                                                   levels = colnames(input_mat_order))
 
 heatmap_annotation_map_selected_order <- selected %>%
   dplyr::select(-ichorcna_tf) %>%
   column_to_rownames(var = "primary")
 
 
 # customize annotation legends
 
 my_annotation_legend_param = list(
   'Age' = list(
     title = "Age",
     at = c(20, 40, 60, 80, 100)
   ),
   'Binary Cohort' = list(
     title = "Binary Cohort",
     at = levels(heatmap_annotation_map$`Binary Cohort`),
     legend_gp = gpar(fill = c("green", "red"))
     
   )
   
   
 )
 
 if(!show_detailed_sample_annotation){
   heatmap_annotation_map_selected_order <- NULL
   my_annotation_legend_param <- NULL
 }
 
 anno <- HeatmapAnnotation(
                          #  tMAD = mad_anno,
                           TF = tf_anno,
                           df = heatmap_annotation_map_selected_order,
                           show_annotation_name = show_annotation_name,
                           annotation_name_side = annotation_name_side,
                           gap = unit(0.15, "cm"),
                           simple_anno_size =  unit(0.3, "cm"),
                           annotation_legend_param = my_annotation_legend_param,
                           border = FALSE
                               )
 
 p <- Heatmap(input_mat_order, 
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              row_title = row_title,
              column_title = col_title,
              # row_title_gp = gpar(fontsize = 6),
              # column_title_gp = gpar(fontsize = 7),
              show_heatmap_legend = show_heatmap_legend,
              heatmap_legend_param = list(title = legend_title, 
                                          at = legend_at
                                          # legend_height = unit(2, "cm"),
                                          # grid_width = unit(0.13, "cm"),
                                          # title_gp = gpar(fontsize = 5),
                                          # labels_gp = gpar(fontsize = 5)
                                          ),
              use_raster = FALSE,
              raster_device = c("png"),
              raster_quality = 30,
              left_annotation = chr_name,
              top_annotation = anno
 )
  
  # add tmad as point distribution (session 3.6 in complexheatmap)
  return(p)
  
}


# segmentation function



seg <- function(x){
  #bicohort <- x$bicohort
  #rowname <- x$rowname
  #primary <- x$primary
  #log2ratio <- x$value
  
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
  
  #ans <- tibble(primary = primary, 
  #              bicohort = bicohort, 
  #              rowname = rowname,
  #              seg = seg_value,
  #              log2ratio = log2ratio )
  
  ans <- x %>% mutate(
    value = seg_value
  )
  
  return(ans)
  
  
}


tMAD_func <- function(x){ 
  abs(mad(x = x, center=0, na.rm = TRUE))
}
