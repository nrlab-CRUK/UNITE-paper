##############################################################################
# parameters
##############################################################################

input_file <- "/Users/wang04/Documents/libs/github/ulyses/data/finaleDB/EE88244.hg19.frag.tsv.bgz.GRanges.rds.1M.rds"

input_file <- "/Users/wang04/Documents/libs/github/ulyses/data/delfi/DL-00001.DL018-DL556.HiSeq_2000_2500.bam.0.1x.mrkdup.bam"

###############################################################################

bam_id <- basename(input_file)

if(stringr::str_detect(input_file, "EE\\d+\\.hg")){
  bam_id <- stringr::str_extract(bam_id, "(^EE\\d+)\\.hg\\d\\d", group = 1)
}

image_file <- paste(input_file, ".ulyses_image.rds", sep = '')
theme_file <- "~/Documents/libs/github/ulyses/explore_packaged/functions/themes.r"

specimen <- "plasma"
which_channel <- "n_isize"
which_bin_size <- "n_100kb_bins_50" 

isize_filter_min <- 90 
isize_filter_max <- 220 


short_range <- if(specimen == "plasma") {
  as.vector(100:150)
} else {as.vector(50:90)}

long_range <- if(specimen == "plasma") { 
  as.vector( 151:220)
} else {as.vector(130:170)}

short <- paste("isize_", short_range, sep = "")
short_label <- paste(min(short_range), "-", max(short_range), sep = '')
long <- paste("isize_", long_range, sep = "")
long_label <- paste( min(long_range), "-", max(long_range), sep = '')
  

ulyses_plot_filename <- paste(input_file, which_channel, which_bin_size, "ulyses_plot.pdf", sep = ".")
ulyses_grid_filename <- paste(input_file, which_channel, which_bin_size, "ulyses_grid.pdf", sep = ".")
ulyses_obj_filename <- paste(input_file, which_bin_size, "ulyses_obj.rds", sep = ".")

##############################################################################
# libs
##############################################################################

library(tidyverse)
library(patchwork)
library(EBImage)
library(gridExtra)
library(grid)
library(ggplotify)
library(reshape2)
library(plotly)



plot_channel <- function(x) {
  
  
  x <- EBImage::normalize(x) %>% reshape2::melt() 
  
  result <-ggplot2::ggplot(x, aes(Var2, Var1)) +
    geom_raster(aes(fill=value), interpolate = FALSE) +
    #scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          legend.key.size = unit(0.3, "cm"),
          legend.text = element_text(size=5),
          legend.position="right",
          legend.title = element_blank())
  
  return(result)
}



# conver long tibble summary to image array
long_to_array <- function(obj) {
  
  channel_names <- unique(obj$channel)
  row_names <- unique(obj$id)
  col_names <- unique(obj$frag_len) %>% 
    paste("isize_", ., sep = "") 
  
  target_dim_names <- list(row_names, col_names, channel_names)
  
  target_dim <- c(length(row_names), 
                  length(col_names), 
                  length(channel_names))
  
  
  image_array <- base::array(obj$pixel, 
                             dim = target_dim, 
                             dimnames = target_dim_names)
  
  return(image_array)
}


##############################################################################
# pre processing
###############################################################################
source(theme_file)
image_array_all <- readRDS(image_file)

image_array <- image_array_all[[which_bin_size]]
image_array_original_isize_range <- dimnames(image_array)[[2]] %>% 
  stringr::str_remove("isize_") %>%
  as.numeric()


if(isize_filter_min < min(image_array_original_isize_range)){
  
  isize_filter_min <- min(image_array_original_isize_range)
  message("'isize_filter_min' is out of range, replaced with range min:")
  message(min(image_array_original_isize_range))
}

if(isize_filter_max > max(image_array_original_isize_range)){
  isize_filter_max <- max(image_array_original_isize_range)
  message("'isize_filter_max' is out of range, replaced with range max.")
  message(max(image_array_original_isize_range))
}

isize_axis_breaks <- c(isize_filter_min, 100, 166, 200, isize_filter_max)
isize_axis_labels <- as.character(isize_axis_breaks)
isize_filter <- paste("isize_", isize_filter_min:isize_filter_max, sep = "")


image_array <- image_array[,isize_filter, ]

bin_start_index <- 1
bin_end_index <- dim(image_array)[1]
bin_axis_breaks <- as.vector(c(bin_start_index, bin_end_index))
bin_axis_labels <- as.character(bin_axis_breaks) 

bin_GC <- apply(image_array[, , "bin_mean_GC"], 1, unique)
bin_mappability <- apply(image_array[, , "bin_mean_mappability"], 1, unique)


###############################################################################
#  channel normalization 
###############################################################################

# motif channels


# NM channels
n_bin <- dimnames(image_array)[[1]] %>% length()
isize_vector <-  isize_filter_min : isize_filter_max
length_array <- replicate(n_bin, isize_vector) %>% 
  t()

#n_bases_array <- image_array[, , "n_isize"] * as.vector(length_array)
## extract nm layers


# pos nm
n_pos_nm_channel <- image_array[, , "n_pos_nm"]
n_pos_nm_channel_per_base <- n_pos_nm_channel / length_array
n_pos_nm_channel_per_base_norm <- n_pos_nm_channel_per_base %>% EBImage::normalize()

# neg nm
n_neg_nm_channel <- image_array[, , "n_neg_nm"]
n_neg_nm_channel_per_base <- n_neg_nm_channel / length_array
n_neg_nm_channel_per_base_norm <- n_neg_nm_channel_per_base %>% EBImage::normalize()


tmp <- EBImage::abind(n_pos_nm_channel_per_base_norm, n_pos_nm_channel %>% EBImage::normalize(), along = 1)

display(tmp)


###############################################################################
# ulyses C/A ratio 
###############################################################################


s1_C_bin_sum <- apply(image_array[, , "n_motif_s1_C"], 1, sum)
s1_G_bin_sum <- apply(image_array[, , "n_motif_s1_G"], 1, sum)
s1_A_bin_sum <- apply(image_array[, , "n_motif_s1_A"], 1, sum)
s1_T_bin_sum <- apply(image_array[, , "n_motif_s1_T"], 1, sum)


s1_C_bin_sum_l_model <- loess(s1_C_bin_sum ~ bin_GC * bin_mappability)
s1_G_bin_sum_l_model <- loess(s1_G_bin_sum ~ bin_GC * bin_mappability)
s1_A_bin_sum_l_model <- loess(s1_A_bin_sum ~ bin_GC * bin_mappability)
s1_T_bin_sum_l_model <- loess(s1_T_bin_sum ~ bin_GC * bin_mappability)

s1_C_bin_sum_corrected <- s1_C_bin_sum - s1_C_bin_sum_l_model$fitted + median(s1_C_bin_sum)
s1_G_bin_sum_corrected <- s1_G_bin_sum - s1_G_bin_sum_l_model$fitted + median(s1_G_bin_sum)
s1_A_bin_sum_corrected <- s1_A_bin_sum - s1_A_bin_sum_l_model$fitted + median(s1_A_bin_sum)
s1_T_bin_sum_corrected <- s1_T_bin_sum - s1_T_bin_sum_l_model$fitted + median(s1_T_bin_sum)

s1_CA_ratio <- s1_C_bin_sum / s1_A_bin_sum
s1_CT_ratio <- s1_C_bin_sum / s1_T_bin_sum
s1_CA_ratio_corrected <- s1_C_bin_sum_corrected / s1_A_bin_sum_corrected
s1_CT_ratio_corrected <- s1_C_bin_sum_corrected / s1_T_bin_sum_corrected



###############################################################################
# ulyses grid plot
###############################################################################
n_channels <- dim(image_array)[3]

channel_names <- as.list(1:n_channels)  %>% setNames(dimnames(image_array)[[3]])

channel_list <- purrr::map(channel_names, .f = EBImage::getFrame, y = image_array )

plot_list <- purrr::map(channel_list, plot_channel)

for (i in 1:n_channels) {
  plot_list[[i]] <- plot_list[[i]] + labs(title = names(plot_list[i]))
}

plot_grid <- purrr::reduce(.x = plot_list, .f = `+`)

ggsave(filename = ulyses_grid_filename, plot = plot_grid, width = 12, height = 10)
message(paste("saved", ulyses_grid_filename))
yulab.utils::o(ulyses_grid_filename)



###############################################################################
# extract summary metrics
###############################################################################

# target_layer matrix

target_matrix_raw <- image_array[, , which_channel]

target_matrix_frac <- target_matrix_raw/sum(target_matrix_raw) 

target_matrix_frac_norm <- target_matrix_frac %>% 
  EBImage::normalize()




# col sum

col_sum <- apply(target_matrix_frac, 2, sum)

#  sd
col_sd <- apply(target_matrix_frac, 2, sd)

col_sd_df <- tibble::enframe(col_sd, name = "isize", value = "col_sd") %>%
  dplyr::mutate(col_sum = col_sum ) %>%
  dplyr::mutate(isize_label = as.numeric(gsub("isize_", "", isize)))

isize_l <- loess(col_sd ~ isize_label, data = col_sd_df, span = 0.2)

col_sd_df_final <- col_sd_df %>%
  dplyr::mutate(col_sd_detrend = col_sd - predict(isize_l))


# cor between bin_GC and count for each isize

col_cor  <- apply(target_matrix_raw, 2, stats::cor, y = bin_GC)


# bin CNV z score
bin_sum <- apply(target_matrix_raw, 1, sum)
l_model <- loess(bin_sum ~ bin_GC * bin_mappability, span = 0.75)
l_model_gc <- loess(bin_sum ~ bin_GC , span = 0.75)
l_model_mappability <- loess(bin_sum ~ bin_mappability, span = 0.85)

bin_sum_gc_fit <- predict(l_model_gc)
bin_sum_mappability_fit <- predict(l_model_mappability)
bin_sum_fit <- predict(l_model_gc)
bin_sum_corrected <- bin_sum -bin_sum_fit + median(bin_sum)
zscore  <- (bin_sum_corrected - mean(bin_sum_corrected))/sd(bin_sum_corrected)


# short 
short_sum <- apply(image_array[, short, which_channel], 1, sum )
short_l_model <- loess(short_sum ~ bin_GC * bin_mappability)
short_sum_corrected <- short_sum - short_l_model$fitted + median(short_sum)

# long 
long_sum <- apply(image_array[, long, which_channel], 1, sum )
long_l_model <- loess(long_sum ~ bin_GC * bin_mappability)
long_sum_corrected <- long_sum - long_l_model$fitted + median(long_sum)

# ratio
sl_ratio <- short_sum/long_sum
sl_ratio_corrected <- short_sum_corrected/long_sum_corrected


bin_df_final <- tibble(bin = names(bin_sum), 
             bin_GC = bin_GC,
             bin_mappability = bin_mappability,
             bin_sum = bin_sum, 
             bin_sum_fit = bin_sum_fit,
             bin_sum_gc_fit = bin_sum_gc_fit,
             bin_sum_mappability_fit = bin_sum_mappability_fit,
             bin_sum_corrected = bin_sum_corrected,
             zscore_corrected = zscore,
             short_sum= short_sum,
             short_sum_corrected = short_sum_corrected,
             long_sum = long_sum,
             long_sum_corrected = long_sum_corrected,
             sl_ratio = sl_ratio,
             sl_ratio_corrected = sl_ratio_corrected,
             s1_CA_ratio = s1_CA_ratio,
             s1_CT_ratio = s1_CT_ratio,
             s1_CA_ratio_corrected = s1_CA_ratio_corrected,
             s1_CT_ratio_corrected = s1_CT_ratio_corrected
             ) %>%
  dplyr::mutate(bin_index = row_number())

bin_df_final$bin <- factor(bin_df_final$bin, levels = names(bin_sum))




###############################################################################
# viz 
###############################################################################


# channels

#im <- as.Image(ulyses_obj$image) %>% 
#  normalize() %>%
#  rotate(-90)

#img_resize = resize(im, w=256, h=256)

#img_t = transpose(img_resize)

#img_neg = max(im) - im

#display( img_neg, method = "raster", bg = "white", frame = 1)

#display( im, method = "raster", all = TRUE)

#writeImage(img_neg, "sample.jpeg", quality = 100)

# channels plus marginal stats 



target_matrix_melt <-target_matrix_frac_norm %>%
  reshape2::melt() 
  #tibble::as_tibble() %>%
  #dplyr::filter(Var2 %in% paste("isize_", isize_filter_min:isize_filter_max, sep = ''))


c1 <- ggplot(target_matrix_melt, aes(Var2, Var1)) +
  geom_raster(aes(fill=value), interpolate = FALSE) +
  #scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size=5),
        legend.position="right",
        legend.title = element_blank())

# count corrected
c2 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, bin_sum_corrected), color = "grey", size = 0.3) +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  coord_flip()+
  theme_row_mid()


# zscore

c3 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, zscore ), color = "orange3", size = 0.3) +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  coord_flip()+
  theme_row_mid()

# n short corrected

c4 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, short_sum_corrected ), color = "grey", size = 0.3) +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  coord_flip() +
  theme_row_mid()




# n long corrected

c5 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, long_sum_corrected ), color = "grey3", size = 0.3) +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  coord_flip() +
  theme_row_mid()

# C/A ratio

c6.0 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, s1_CT_ratio_corrected ), color = "blue3", size = 0.3) +
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = bin_axis_breaks,
                     labels = bin_axis_labels,
                     expand = c(0, 0)) + 
  labs(x = paste("Bin", " (", which_bin_size, ") "), y = paste("C/T Ratio Corrected")) +
  coord_flip() +
  theme_row_mid()




# s/l ratio from corrected s and l count

c6 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, sl_ratio_corrected ), color = "blue3", size = 0.3) +
  scale_y_continuous(position = "right") +
  scale_x_continuous(breaks = bin_axis_breaks,
                     labels = bin_axis_labels,
                     expand = c(0, 0)) + 
  labs(x = paste("Bin", " (", which_bin_size, ") "), y = paste("S/L", "\n",short_label, "/", long_label, " bp", sep = "")) +
  coord_flip() +
  theme_row_border()


# raw fragmentation profile

c7 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sum, group = 1), size = 0.7, color = "grey") +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.05), expand = c(0, 0), position = "right") +
  labs(y = "Fraction of count") +
  theme_col_mid()


# raw sd line
c8 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sd, group = 1), size = 0.7, color = "orange3" ) +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.000015), expand = c(0, 0), position = "right") +
  #geom_hline(yintercept = max(col_sd_df_final$col_sd)) +
  labs(y = "SD Scaled") +
  theme_col_mid()
  

# detrended  sdline
c9 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sd_detrend, group = 1), color = "grey", size = 0.7) +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), 
                     breaks = isize_axis_breaks, 
                     labels = isize_axis_labels, 
                     expand = c(0, 0)) + 
  scale_y_continuous( position = "right") +
  labs(x = "Fragment length (bp)", y = "SD Scaled Detrended") +
  theme_col_border()


ulyses_plot <- c6 + c6.0 + c3 + c2 + c1 + 
  plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + c7 +
  plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + c8 +
  plot_spacer() + plot_spacer() + plot_spacer() +plot_spacer() + c9 +
  plot_layout(ncol = 5, nrow = 4, widths = c(1, 1, 1, 1, 6), heights = c(6, 1, 1, 1))



ulyses_plot


ggsave(filename = ulyses_plot_filename, plot = ulyses_plot, width = 10, height = 8 )
message(paste("Saved", ulyses_plot_filename))
yulab.utils::o(ulyses_plot_filename)





###############################################################################
# QC plot of CNV and zscore
###############################################################################




# GC-count
p1 <- ggplot(bin_df_final) +
  geom_point(aes(bin_GC, bin_sum), color = "black", size = 0.3) +
  geom_line(aes(bin_GC, bin_sum_gc_fit), color = "orange3") +
  theme_classic()

p2 <- ggplot(bin_df_final) +
  geom_point(aes(bin_GC, bin_sum_corrected), color = "grey", size = 0.3) +
  theme_classic()

pc1 <- (p1 / p2)


p_ranges_x <- c(ggplot_build(pc1[[1]])$layout$panel_scales_x[[1]]$range$range,
                ggplot_build(pc1[[2]])$layout$panel_scales_x[[1]]$range$range)

p_ranges_y <- c(ggplot_build(pc1[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(pc1[[2]])$layout$panel_scales_y[[1]]$range$range)

pc1_final <- pc1 & 
  xlim(min(p_ranges_x), max(p_ranges_x)) & 
  ylim(min(p_ranges_y), max(p_ranges_y))



# QC viz
# mappability-count
p3 <- ggplot(bin_df_final) +
  geom_point(aes(bin_mappability, bin_sum), size = 0.3) +
  geom_line(aes(bin_mappability, bin_sum_mappability_fit), color = "orange3") +
  theme_classic()

p4 <- ggplot(bin_df_final) +
  geom_point(aes(bin_mappability, bin_sum_corrected), color = "grey", size = 0.3) +
  theme_classic()

pc2 <- (p3 / p4)


p_ranges_x <- c(ggplot_build(pc2[[1]])$layout$panel_scales_x[[1]]$range$range,
                ggplot_build(pc2[[2]])$layout$panel_scales_x[[1]]$range$range)

p_ranges_y <- c(ggplot_build(pc2[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(pc2[[2]])$layout$panel_scales_y[[1]]$range$range)

pc2_final <- pc2 & 
  xlim(min(p_ranges_x), max(p_ranges_x)) & 
  ylim(min(p_ranges_y), max(p_ranges_y))


# count of bin

p5 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, bin_sum), color = "black", size = 0.3) +
  theme_classic() 

p6 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, bin_sum_corrected), color = "grey", size = 0.3) +
  theme_classic() 



pc3 <- (p5 / p6)

p_ranges_x <- c(ggplot_build(pc3[[1]])$layout$panel_scales_x[[1]]$range$range,
                ggplot_build(pc3[[2]])$layout$panel_scales_x[[1]]$range$range)

p_ranges_y <- c(ggplot_build(pc3[[1]])$layout$panel_scales_y[[1]]$range$range,
                ggplot_build(pc3[[2]])$layout$panel_scales_y[[1]]$range$range)

pc3_final <- pc3 & 
  xlim(min(p_ranges_x), max(p_ranges_x)) & 
  ylim(min(p_ranges_y), max(p_ranges_y))

p7 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, zscore ), color = "orange3", size = 0.3) +
  theme_classic() 

p_final <- (pc1_final | pc2_final) / (pc3_final /p7)

ggsave(paste0(ulyses_plot_filename, "_CNV_normalization_plot.pdf"), plot = p_final, height = 10, width = 7)

yulab.utils::o(paste0(ulyses_plot_filename, "_CNV_normalization_plot.pdf"))






###############################################################################
# integrate as list 
###############################################################################

ulyses_obj <- list("image" = image_array,
                   "col_sd" = col_sd,
                   "short_sum_corrected" = short_sum_corrected,
                   "long_sum_corrected" = long_sum_corrected, 
                   "sl_ratio_corrected" = sl_ratio_corrected,
                   "s1_CT_ratio_corrected" = s1_CT_ratio_corrected,
                   "zscore" = bin_df_final$zscore_corrected)

saveRDS(object = ulyses_obj, file = ulyses_obj_filename )


