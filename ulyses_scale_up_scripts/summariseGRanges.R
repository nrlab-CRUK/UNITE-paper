
library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)

args = commandArgs(trailingOnly=TRUE)

inputfile  <- as.character(args[1])

output_dir <- as.character(args[2])

overwrite <-  as.logical(args[3]) 


frag <- readRDS(inputfile)
n_frag <- length(frag)
result <- tibble("file" = inputfile, "n_frag" = n_frag)
filename <- file.path(output_dir, paste0(basename(inputfile), ".summary.csv"))


if( (!file.exists(filename))|| isTRUE(overwrite)){
  write_csv(x = result, file = filename )
  message("csv file written successfully.")
} else {
  message("output file exists. \n")
  message(filename, " skipped. \n")
}


