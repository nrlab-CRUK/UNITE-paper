
# Library -----------------------------------------------------------------

suppressMessages(library(tidyverse))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(EBImage))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(plotly))
# library(ggplotify)
# library(plotly)
#healthy: EE87922.hg19.frag.tsv.bgz.GRanges.rds.1M.rds.ulyses_image.rds
#cancer: EE87033.hg19.frag.tsv.bgz.GRanges.rds.1M.rds.ulyses_image.rds
library(argparser)

# Define command line arguments
parser <- arg_parser(description = "Prepare tensor and plots.")
parser <- add_argument(parser, "--image_file", type = "character", 
                       default = "/Users/wang04/ulyses/5_ulyses_plot/EE87922.hg19.frag.tsv.bgz.GRanges.rds.1M.rds.ulyses_image.rds",
                       #default = "/Users/wang04/ulyses/5_ulyses_plot/EE87033.hg19.frag.tsv.bgz.GRanges.rds.1M.rds.ulyses_image.rds",
                       help = "Path to the image file.")
parser <- add_argument(parser, "--specimen", type = "character", 
                       help = "Specimen name.", default = "plasma")
parser <- add_argument(parser, "--which_channel", type = "character", 
                       help = "Channel name.", default = "n_isize")
parser <- add_argument(parser, "--which_bin_size", type = "character", 
                       help = "Bin size name.", default = "n_100kb_bins_50")
parser <- add_argument(parser, "--isize_filter_min", type = "numeric", 
                       help = "Minimum isize filter.", default = 100)
parser <- add_argument(parser, "--isize_filter_max", type = "numeric", 
                       help = "Maximum isize filter.", default = 220)
parser <- add_argument(parser, "--theme_file", type = "character", 
                       help = "Theme file.", 
                       default = "~/ulyses/explore_packaged/functions/themes.r")

# add a param called ciridis_color_name
parser <- add_argument(parser, "--viridis_color_name", type = "character", 
                       help = "Viridis color name.", default = "turbo")
parser <- add_argument(parser, "--viridis_direction", type = "numeric",
                       help = "Viridis direction.", default = 1)
parser <- add_argument(parser, "--axis_linewidth", type = "numeric",
                       help = "Axis linewidth.", default = 0.15)
parser <- add_argument(parser, "--viridis_limits", type = "numeric", nargs = 2,
                       help = "Viridis limits.", default = c(-4, 4))
parser <- add_argument(parser, "--viridis_limits_breaks", type = "numeric", nargs = 5,
                       help = "Viridis limits breaks.", default = c(-4, -2, 0, 2, 4))


# Parse command line arguments
args <- parse_args(parser)

# Accessing values from parsed arguments
image_file <- args$image_file
specimen <- args$specimen
which_channel <- args$which_channel
which_bin_size <- args$which_bin_size
isize_filter_min <- args$isize_filter_min
isize_filter_max <- args$isize_filter_max
theme_file <- args$theme_file


viridis_color_name <- args$viridis_color_name
viridis_direction <- args$viridis_direction
axis_linewidth <- args$axis_linewidth
viridis_limits <- args$viridis_limits
viridis_limits_breaks <- args$viridis_limits_breaks


# preprocessing 
source(theme_file)

# get bam_id , useful for meta file
if(stringr::str_detect(image_file, "EE\\d+\\.hg")){
  bam_id <- stringr::str_extract(basename(image_file), "(^EE\\d+)\\.hg\\d\\d", group = 1)
} else {
  bam_id <- stringr::str_extract(basename(image_file), "(^(SLX|DL)\\-\\d+\\.\\w+\\d+(\\w+)?(\\-\\w+\\d+(\\w+)?)?)\\.", group = 1)
}

short_range <- if(specimen == "plasma") {
  as.vector(100:150)
} else {as.vector(50:90)}

long_range <- if(specimen == "plasma") { 
  as.vector( 151:220)
} else {as.vector(130:170)}

# prepare short long labels in the Ulyses plot
short <- paste("isize_", short_range, sep = "")
short_label <- paste(min(short_range), "-", max(short_range), sep = '')
long <- paste("isize_", long_range, sep = "")
long_label <- paste( min(long_range), "-", max(long_range), sep = '')

# output file names
ulyses_plot_filename <- paste(image_file, which_channel, which_bin_size, "ulyses_plot.pdf", sep = ".")
ulyses_grid_filename <- paste(image_file, which_channel, which_bin_size, "ulyses_grid.pdf", sep = ".")
ulyses_obj_filename <- paste(image_file, which_bin_size, "ulyses_obj.rds", sep = ".")



# Functions ---------------------------------------------------------------


plot_channel_archived <- function(x, 
                         normalize = FALSE) {
  
  if(normalize) {
    x <- EBImage::normalize(x) 
  }
  
  x <- reshape2::melt(x, varnames = c("bin", "isize")) 
  
  result <-ggplot2::ggplot(x, aes(isize, bin)) +
    geom_raster(aes(fill=value), interpolate = FALSE) +
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

plot_channel <- function(x, 
                         normalize = FALSE) {
  
  if(normalize) {
    x <- EBImage::normalize(x) 
  }
  
  x <- reshape2::melt(x, varnames = c("bin", "isize")) 
  
  result <- ggplot2::ggplot(x, aes(isize, bin, fill = value)) +
    geom_tile() +
    scale_fill_gradient(limits = c(0, 1)) +
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


# Convert long tibble summary to image array
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

matrix_frac_norm <- function(x, 
                             which_channel, 
                             frac = TRUE, 
                             normalize = TRUE) {
  ans <- x[, , which_channel]
  
  if(frac) {
    ans <- target_matrix_raw/sum(target_matrix_raw) 
  }
  
  if(normalize) {
    ans <- EBImage::normalize(ans)
  }
  
  return(ans)
  
}


motif_total_frac_norm <- function(x, 
                                  which_channel, 
                                  normalize = FALSE) {
  n_total <- sum(x[,,"n_isize"]) 
  
  target_matrix_raw <- x[, , which_channel]
  
  ans <- target_matrix_raw/n_total 
  
  if(normalize) {
    ans <- EBImage::normalize(ans)
  }
  
  return(ans)
  
}


nm_norm_func <- function(array, 
                         per_base = TRUE,
                         normalize = TRUE
){
  if(per_base) {
    ans <- array / length_array
  }
  
  if(normalize) {
    ans <- ans %>% EBImage::normalize()
  }
}

# Pre-processing ----------------------------------------------------------


image_array_all <- readRDS(image_file)

image_array <- image_array_all[[which_bin_size]]

# remove problematic bins called "6p_6"
which_bin <- which(dimnames(image_array)[[1]] == "6p_6")
image_array <- image_array[-which_bin,,]



image_array_original_isize_range <- dimnames(image_array)[[2]] %>% 
  # remove element containing "NA"
  purrr::keep(.p = ~!stringr::str_detect(.x, "NA")) %>%
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

# bin GC and mappability

bin_GC <- apply(image_array[, , "bin_mean_GC"], 1, mean, na.rm = TRUE)
bin_mappability <- apply(image_array[, , "bin_mean_mappability"], 1, mean, na.rm = TRUE)

# make an array, n cols equal to dim(image_array[,,"bin_mean_GC"])[2], each col is bin_GC
bin_GC_array <- replicate(dim(image_array[,,"bin_mean_GC"])[2], bin_GC)
bin_mappability_array <- replicate(dim(image_array[,,"bin_mean_mappability"])[2], bin_mappability)

# substitute bin_GC and bin_mappability with array
image_array[,,"bin_mean_GC"] <- bin_GC_array
image_array[,,"bin_mean_mappability"] <- bin_mappability_array



# Channel normalization ---------------------------------------------------


# normalize n_isize and n_motif layers together
channel_to_keep <- dimnames(image_array)[[3]] %>% 
  purrr::keep(.p = ~stringr::str_detect(.x,"n_isize|n_motif"))
image_array_keep <- image_array[,, channel_to_keep]
# the seperate = FALSE is important, otherwise the normalization will be done per channel
image_array_keep_norm <- EBImage::normalize(image_array_keep, separate = FALSE)
# array to list
isize_motif_channel_list <- purrr::array_branch(image_array_keep_norm, margin = 3) 



# Normalize NM channel
nm_channel <- dimnames(image_array)[[3]] %>% 
  purrr::keep(.p = ~stringr::str_detect(.x, "nm"))

if(length(nm_channel) > 0) {
  
  isize_vector <-   dimnames(image_array) %>%
    purrr::pluck(2) %>%
    stringr::str_remove(pattern = "^isize_") %>%
    as.numeric()
  
  length_array <- replicate(dim(image_array)[[1]], isize_vector) %>% t()
  
  nm_channel_array <- image_array[,, nm_channel]
  nm_channel_list <- purrr::array_branch(nm_channel_array, margin = 3) %>%
    purrr::map(.f = nm_norm_func, per_base = TRUE, normalize = TRUE )
  
  # nm_channel_list <- purrr::map(seq(1, length(nm_channel), 1), 
  #                               .f = EBImage::getFrame, 
  #                               y = nm_channel_array) %>%
  #   setNames(nm_channel) %>%
  #   purrr::map(.f = nm_norm_func, per_base = TRUE, normalize = TRUE )
  
} else {
  nm_channel_list = NULL
  message("NM layers do not exist, skipping...")
}



# normalize housekeeping layers
housekeeping_channel <- dimnames(image_array)[[3]] %>% 
  purrr::keep(.p = ~stringr::str_detect(.x,"^bin"))

# if there is no housekeeping channel, skip
if(length(housekeeping_channel) > 0) {
  housekeeping_channel_array <- image_array[, , housekeeping_channel]
  
  housekeeping_channel_list <- purrr::array_branch(housekeeping_channel_array, margin = 3) %>%
    purrr::map(.f = EBImage::normalize)
  
} else {
  housekeeping_channel_list = NULL
  message("Housekeeping layers do not exist, skipping...")
}




# combine all 
channel_list_combined <- c( isize_motif_channel_list, 
                            nm_channel_list, 
                            housekeeping_channel_list) 

tensor <- EBImage::combine(channel_list_combined) 

dimnames(tensor)[[3]] <- names(channel_list_combined)




# n_isize channel feature extract -----------------------------------------


target_matrix_raw <- matrix_frac_norm(image_array, 
                                      which_channel = "n_isize",
                                      frac = FALSE,
                                      normalize = FALSE)

target_matrix_frac <- matrix_frac_norm(image_array, 
                                       which_channel = "n_isize",
                                       frac = TRUE,
                                       normalize = FALSE)

target_matrix_frac_norm <- matrix_frac_norm(image_array, 
                                            which_channel = "n_isize",
                                            frac = TRUE,
                                            normalize = TRUE)

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
cnv_log2Ratio <- log2(bin_sum_corrected/median(bin_sum_corrected))

# short 
short_sum <- apply(image_array[, short, which_channel], 1, sum )
short_l_model <- loess(short_sum ~ bin_GC * bin_mappability)
short_sum_corrected <- short_sum - short_l_model$fitted + median(short_sum)

# long 
long_sum <- apply(image_array[, long, which_channel], 1, sum )
long_l_model <- loess(long_sum ~ bin_GC * bin_mappability)
long_sum_corrected <- long_sum - long_l_model$fitted + median(long_sum)

# s/l ratio
sl_ratio <- short_sum/long_sum
sl_ratio_corrected <- short_sum_corrected/long_sum_corrected


# motif ratio

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


# gather the results
bin_df_final <- tibble(bin = names(bin_sum), 
                       bin_GC = bin_GC,
                       bin_mappability = bin_mappability,
                       bin_sum = bin_sum, 
                       bin_sum_fit = bin_sum_fit,
                       bin_sum_gc_fit = bin_sum_gc_fit,
                       bin_sum_mappability_fit = bin_sum_mappability_fit,
                       bin_sum_corrected = bin_sum_corrected,
                       zscore_corrected = zscore,
                       cnv_log2Ratio = cnv_log2Ratio,
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


# Ulyses n_isize plot -----------------------------------------------------

# channels
target_matrix_melt <-target_matrix_frac_norm %>%
  reshape2::melt(varnames = c("bin", "isize")) 



# count corrected
c2 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, bin_sum_corrected), color = "grey", size = 0.3) +
  scale_y_continuous(position = "right", expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  coord_flip()+
  theme_row_mid()



# zscore

c3 <- ggplot(bin_df_final) +
  geom_point(aes(bin_index, zscore_corrected ), color = "orange3", size = 0.3) +
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



# raw sd line
c8 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sd, group = 1), linewidth = 0.3, color = "orange" ) +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), 
                     breaks = isize_axis_breaks, 
                     labels = isize_axis_labels, 
                     expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.000030), expand = c(0, 0), position = "left") +
  #geom_hline(yintercept = max(col_sd_df_final$col_sd)) +
  labs(y = "SD (Scaled)",
       x ="Fragment Length (bp)") +
  theme_col_mid()



# detrended  sdline
c9 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sd_detrend, group = 1), color = "grey", linewidth = 0.7) +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), 
                     breaks = isize_axis_breaks, 
                     labels = isize_axis_labels, 
                     expand = c(0, 0)) + 
  scale_y_continuous( position = "right") +
  labs(x = "Fragment length (bp)", y = "SD Scaled Detrended") +
  theme_col_border()

################################################################################
# useful panels
################################################################################

c1 <- ggplot(target_matrix_melt, aes(isize, bin)) +
  geom_raster(aes(fill=value), interpolate = FALSE) +
  scale_fill_viridis_c(option =viridis_color_name, direction = viridis_direction) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size=5),
        legend.position="right",
        legend.title = element_blank())

((
  
c1_no_legend <- ggplot(target_matrix_melt, aes(isize, bin)) +
  geom_raster(aes(fill=value), interpolate = FALSE) +
  # set color to viridis plasma
  scale_fill_viridis_c(option = viridis_color_name, direction = viridis_direction) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size=5),
        legend.position="none",
        legend.title = element_blank())
  
))

# s/l ratio from corrected s and l count

input <- bin_df_final %>% mutate(sl_ratio_corrected_z = scale(sl_ratio_corrected))

c6 <- ggplot(input) +
  geom_point(aes(x = bin_index, 
                 y = sl_ratio_corrected, 
                 color = sl_ratio_corrected_z), 
             size = 0.1) +
  scale_y_continuous(position = "left", sec.axis = sec_axis(~ scale(.), name = "z-score")) +
  scale_x_continuous(breaks = bin_axis_breaks,
                     labels = bin_axis_labels,
                     expand = c(0, 0)
                     ) + 
  #labs(x = paste("Genomic Bins", " (", which_bin_size, ") "), 
  #     y = paste("S/L", "\n",short_label, "/", long_label, " bp", sep = "")) +
  labs(x = paste("Genomic Bin Index", " (5Mb) "), 
       y = paste("S/L", "\n",short_label, "/", long_label, " bp", sep = "")) +
  # set y title to 'S/L Ratio'
  ylab("S/L Ratio") +
  coord_flip() +
  scale_color_viridis_c(option = viridis_color_name, direction = viridis_direction, limits= viridis_limits, breaks = viridis_limits_breaks) +
  scale_fill_viridis_c(option = viridis_color_name, direction = viridis_direction, limits= viridis_limits, breaks = viridis_limits_breaks) +
  theme_row_border() +
  # set axis line width = axis_line_width
  theme(axis.line = element_line(linewidth = axis_linewidth)) +
  # set axis tick thickness  = axis_line_width
  theme(axis.ticks = element_line(linewidth = axis_linewidth)) +
  theme(legend.position = "top") +
  # remove legend title
  theme(legend.title = element_blank()) +
  # make legend smaller
  theme(legend.key.size = unit(0.18, "cm")) +
  # legend text size to 5
  theme(legend.text = element_text(size = 5)) +
  # rotate legend text 45 degrees
  theme(legend.text = element_text(angle = 45)) 

# C/T ratio

c6.0 <- ggplot(bin_df_final %>% 
                 mutate(tc_ratio = 1/s1_CT_ratio_corrected) %>%
                 mutate(tc_ratio_z = scale(tc_ratio))) +
  geom_point(aes(bin_index, 
                 tc_ratio , 
                 color = tc_ratio_z), 
             size = 0.1) +
  scale_y_continuous(position = "left", sec.axis = sec_axis(~ scale(.), name = "z-score")) +
  scale_x_continuous(breaks = bin_axis_breaks,
                     labels = bin_axis_labels,
                     expand = c(0, 0)) + 
  labs(x = paste("Bin", " (", which_bin_size, ") "), y = paste("T/C Ratio")) +
  coord_flip() +
scale_color_viridis_c(option = viridis_color_name, direction = viridis_direction, limits= viridis_limits) +
  theme_row_mid() +
  theme(legend.position = "top") +
  # remove legend title
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(linewidth = axis_linewidth)) +
  theme(axis.ticks = element_line(linewidth = axis_linewidth)) +
  # make legend smaller
  theme(legend.key.size = unit(0.18, "cm")) +
  # legend text size to 5
  theme(legend.text = element_text(size = 5)) +
  # rotate legend text 45 degrees
  theme(legend.text = element_text(angle = 45)) 

c_cnv_log2Ratio <- ggplot(bin_df_final %>% 
                            mutate(cnv_log2Ratio_z = scale(cnv_log2Ratio))) +
  geom_point(aes(bin_index, 
                 cnv_log2Ratio, 
                 color = cnv_log2Ratio_z), 
             size = 0.1) +
  scale_y_continuous(position = "left", sec.axis = sec_axis(~ scale(.), name = "z-score")) +
  scale_x_continuous(expand = c(0, 0)) + 
scale_color_viridis_c(option = viridis_color_name, direction = viridis_direction, limits= viridis_limits) +
  # set color manually
  labs(y = "Log2Ratio") +
  coord_flip()+
  theme_row_mid() +
  # legend none
  theme(legend.position = "top") +
  # remove legend title
  theme(legend.title = element_blank()) +
  theme(axis.line = element_line(linewidth = axis_linewidth)) +
  theme(axis.ticks = element_line(linewidth = axis_linewidth)) +
  # make legend smaller
  theme(legend.key.size = unit(0.18, "cm")) +
  # legend text size to 5
  theme(legend.text = element_text(size = 5)) +
  # rotate legend text 45 degrees
  theme(legend.text = element_text(angle = 45)) 

# raw fragmentation profile

c7 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sum, group = 1), linewidth = 0.3, color = "orange") +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), 
                     breaks = isize_axis_breaks, 
                     labels = isize_axis_labels, 
                     expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.05), 
                     breaks = c(0,0.03, 0.05),
                     expand = c(0, 0), 
                     position = "left") +
  labs(y = "Frac.",
       x ="Fragment Length (bp)") +
  theme_col_mid() +
  # no need to rotate y axis label
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5))+
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5))+
  theme(plot.margin = margin(0, 0, 0, 0)) +
  # remove the x axis title, labels and tick completely
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(axis.line = element_line(linewidth = axis_linewidth)) +
  theme(axis.ticks = element_line(linewidth = axis_linewidth)) 


c8_2 <- ggplot(col_sd_df_final) +
  geom_line(aes(isize_label, col_sd, group = 1), linewidth = 0.3, color = "orange" ) +
  scale_x_continuous(limits = c(isize_filter_min, isize_filter_max), 
                     breaks = isize_axis_breaks, 
                     labels = isize_axis_labels, 
                     expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0, 0.000030), expand = c(0, 0), position = "left") +
  #geom_hline(yintercept = max(col_sd_df_final$col_sd)) +
  labs(y = "SD of frac.",
       x ="Fragment Length (bp)") +
  theme_col_border() +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  # y axis titile rotate 90 degree
  theme(axis.title.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
  # reduce the plot margin to 0
  theme(plot.margin = margin(0, 0, 0, 0)) +
  theme(axis.line = element_line(linewidth = axis_linewidth)) +
  theme(axis.ticks = element_line(linewidth = axis_linewidth)) 



  ulyses_n_isize_plot <- c6 + c6.0 + c3 + c_cnv_log2Ratio + c1 + 
    plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + c7 +
    plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() + c8 +
    plot_spacer() + plot_spacer() + plot_spacer() +plot_spacer() + c9 +
    plot_layout(ncol = 5, nrow = 4, widths = c(1, 1, 1, 1, 6), heights = c(6, 1, 1, 1))

((
  
  ulyses_n_isize_plot <- c6 + c6.0 +  c_cnv_log2Ratio +  c1_no_legend + 
    plot_spacer() + plot_spacer() +  plot_spacer()  +  c7 +
    plot_spacer() + plot_spacer() +  plot_spacer() +  c8_2 +
    plot_layout(ncol = 4, nrow = 3, widths = c(1, 1, 1, 5.3), heights = c(5.9, 1, 1))
  
))

ggsave(filename = ulyses_plot_filename, plot = ulyses_n_isize_plot, 
       width = 9, 
       height = 9,
       units = "cm")

########## plot with legend

# plot c1 with legend to an independent file
ulyses_plot_filename_with_legend <- gsub(".pdf", "_with_legend.pdf", ulyses_plot_filename)


((
  
  ulyses_n_isize_plot_with_legend <- c6 + c6.0 +  c_cnv_log2Ratio +  c1 + 
    plot_spacer() + plot_spacer() +  plot_spacer()  +  c7 +
    plot_spacer() + plot_spacer() +  plot_spacer() +  c8_2 +
    plot_layout(ncol = 4, nrow = 3, widths = c(1, 1, 1, 5.3), heights = c(5.9, 1, 1))
  
))

ggsave(filename = ulyses_plot_filename_with_legend, plot = ulyses_n_isize_plot_with_legend, 
       width = 11.5, 
       height = 9,
       units = "cm")

# narrower ---------------------
ulyses_plot_filename_with_legend_narrow <- gsub(".pdf", "_with_legend_narrow.pdf", ulyses_plot_filename)
((
  
  ulyses_n_isize_plot_with_legend_narrow <- c6 + c6.0 +  c_cnv_log2Ratio +  c1 + 
    plot_spacer() + plot_spacer() +  plot_spacer()  +  c7 +
    plot_spacer() + plot_spacer() +  plot_spacer() +  c8_2 +
    plot_layout(ncol = 4, nrow = 3, widths = c(1, 1, 1, 3.2), heights = c(5.9, 1, 1))
  
))

ggsave(filename = ulyses_plot_filename_with_legend_narrow, plot = ulyses_n_isize_plot_with_legend_narrow, 
       width = 10.5, 
       height = 9,
       units = "cm")
message(paste("Saved", ulyses_plot_filename_with_legend_narrow))

# yulabutiles open the pdf file
if (interactive()) {
  utils::browseURL(ulyses_plot_filename_with_legend_narrow)
}

# use ggplotly to save ulyses_n_isize_plot_with_legend_narrow to html
ulyses_plot_filename_with_legend_narrow_html <- gsub(".pdf", ".html", ulyses_plot_filename_with_legend_narrow)
fig1 <- ggplotly(c6) %>% layout(legend = list(orientation = "h"))
fig2 <- ggplotly(c6.0)
fig3 <- ggplotly(c_cnv_log2Ratio)
fig4 <- ggplotly(c1)
fig <- subplot(fig1, fig2, fig3, fig4, widths = c(0.16, 0.16, 0.16, 0.51))
htmlwidgets::saveWidget(fig, ulyses_plot_filename_with_legend_narrow_html, selfcontained = TRUE)
message(paste("Saved", ulyses_plot_filename_with_legend_narrow_html))

# 
# 
# # Ulyses channel grid plot ------------------------------------------------
# 
# plot_list <- purrr::map(channel_list_combined, plot_channel, normalize = FALSE)
# 
# for (i in 1:length(channel_list_combined)) {
#   plot_list[[i]] <- plot_list[[i]] + labs(title = names(plot_list[i]))
# }
# 
#   ulyses_grid_plot <- purrr::reduce(.x = plot_list, .f = `+`)
# 
# ggsave(filename = ulyses_grid_filename, plot = ulyses_grid_plot, width = 12, height = 10)
# message(paste("saved", ulyses_grid_filename))
# 
# # QC plot of CNV  ---------------------------------------------------------
# # GC-count
# p1 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_GC, bin_sum), color = "black", size = 0.3) +
#   geom_line(aes(bin_GC, bin_sum_gc_fit), color = "orange3") +
#   theme_classic()
# 
# p2 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_GC, bin_sum_corrected), color = "grey", size = 0.3) +
#   theme_classic()
# 
# pc1 <- (p1 / p2)
# 
# 
# p_ranges_x <- c(ggplot_build(pc1[[1]])$layout$panel_scales_x[[1]]$range$range,
#                 ggplot_build(pc1[[2]])$layout$panel_scales_x[[1]]$range$range)
# 
# p_ranges_y <- c(ggplot_build(pc1[[1]])$layout$panel_scales_y[[1]]$range$range,
#                 ggplot_build(pc1[[2]])$layout$panel_scales_y[[1]]$range$range)
# 
# pc1_final <- pc1 & 
#   xlim(min(p_ranges_x), max(p_ranges_x)) & 
#   ylim(min(p_ranges_y), max(p_ranges_y))
# 
# 
# # QC viz
# # mappability-count
# p3 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_mappability, bin_sum), size = 0.3) +
#   geom_line(aes(bin_mappability, bin_sum_mappability_fit), color = "orange3") +
#   theme_classic()
# 
# p4 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_mappability, bin_sum_corrected), color = "grey", size = 0.3) +
#   theme_classic()
# 
# pc2 <- (p3 / p4)
# 
# 
# p_ranges_x <- c(ggplot_build(pc2[[1]])$layout$panel_scales_x[[1]]$range$range,
#                 ggplot_build(pc2[[2]])$layout$panel_scales_x[[1]]$range$range)
# 
# p_ranges_y <- c(ggplot_build(pc2[[1]])$layout$panel_scales_y[[1]]$range$range,
#                 ggplot_build(pc2[[2]])$layout$panel_scales_y[[1]]$range$range)
# 
# pc2_final <- pc2 & 
#   xlim(min(p_ranges_x), max(p_ranges_x)) & 
#   ylim(min(p_ranges_y), max(p_ranges_y))
# 
# 
# # count of bin
# 
# p5 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_index, bin_sum), color = "black", size = 0.3) +
#   theme_classic() 
# 
# p6 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_index, bin_sum_corrected), color = "grey", size = 0.3) +
#   theme_classic() 
# 
# 
# 
# pc3 <- (p5 / p6)
# 
# p_ranges_x <- c(ggplot_build(pc3[[1]])$layout$panel_scales_x[[1]]$range$range,
#                 ggplot_build(pc3[[2]])$layout$panel_scales_x[[1]]$range$range)
# 
# p_ranges_y <- c(ggplot_build(pc3[[1]])$layout$panel_scales_y[[1]]$range$range,
#                 ggplot_build(pc3[[2]])$layout$panel_scales_y[[1]]$range$range)
# 
# pc3_final <- pc3 & 
#   xlim(min(p_ranges_x), max(p_ranges_x)) & 
#   ylim(min(p_ranges_y), max(p_ranges_y))
# 
# p7 <- ggplot(bin_df_final) +
#   geom_point(aes(bin_index, zscore_corrected ), color = "orange3", size = 0.3) +
#   theme_classic() 
# 
# ulyses_qc_plot <- (pc1_final | pc2_final) / (pc3_final /p7)
# 
# ggsave(paste0(ulyses_plot_filename, "_CNV_normalization_plot.pdf"), plot = ulyses_qc_plot, height = 10, width = 7)
# 
# 
# # Integrate as a list -----------------------------------------------------
# 
# # tensor filter
# dm <- dimnames(tensor)
# 
# # if "n_neg_nm" and "n_pos_nm" exist in dm
# if("n_neg_nm" %in% dm[[3]] & "n_pos_nm" %in% dm[[3]]) {
#   layers_to_keep <- c(
#   "n_isize",
#   "n_motif_umono3_C",
#   "n_motif_umono3_T",
#   "n_motif_umono2_C",
#   "n_motif_umono2_T",
#   "n_motif_umono1_C",
#   "n_motif_umono1_T",
#   "n_motif_smono1_C",
#   "n_motif_smono1_T",
#   "n_motif_smono2_C",
#   "n_motif_smono2_T",
#   "n_motif_smono3_C",
#   "n_motif_smono3_T",
#   "n_neg_nm",
#   "n_pos_nm",
#   "bin_mean_GC",
#   "bin_mean_mappability"
#   )
# 
# } else {
#   layers_to_keep <- c(
#   "n_isize",
#   "n_motif_umono3_C",
#   "n_motif_umono3_T",
#   "n_motif_umono2_C",
#   "n_motif_umono2_T",
#   "n_motif_umono1_C",
#   "n_motif_umono1_T",
#   "n_motif_smono1_C",
#   "n_motif_smono1_T",
#   "n_motif_smono2_C",
#   "n_motif_smono2_T",
#   "n_motif_smono3_C",
#   "n_motif_smono3_T",
#   "bin_mean_GC",
#   "bin_mean_mappability"
#   )
# 
# }
# 
# 
# 
# tensor <- tensor[,, layers_to_keep]
# 
# ulyses_obj <- list("tensor" = tensor,
#                    "image_file" = image_file,
#                    "bam_id" = bam_id,
#                    "specimen" = specimen,
#                    "which_bin_size" = which_bin_size ,
#                    "isize_filter_min" = isize_filter_min ,
#                    "isize_filter_max" = isize_filter_max ,
#                    "bin_df_final" = bin_df_final,
#                    "col_sd_df_final" = col_sd_df_final,
#                    "ulyses_grid_plot" = ulyses_grid_plot,
#                    "ulyses_n_isize_plot" = ulyses_n_isize_plot,
#                    "ulyses_qc_plot" = ulyses_qc_plot
# )
# 
# saveRDS(object = ulyses_obj, file = ulyses_obj_filename )
# 
# message("saved ulyses obj file: \n", ulyses_obj_filename )





