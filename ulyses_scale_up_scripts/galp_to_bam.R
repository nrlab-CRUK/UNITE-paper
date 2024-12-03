
library(tidyverse)
library(magrittr)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)

args = commandArgs(trailingOnly=TRUE)

inputfile  <- as.character(args[1])

output_dir <- as.character(args[2])

overwrite <-  as.logical(args[3]) 


frag <- readRDS(inputfile)
filename <- file.path(output_dir, paste0(basename(inputfile), ".bam"))
filename_fh <- Rsamtools::BamFile(file = filename, index = filename)

# export frag as bam file
if( (!file.exists(filename))|| isTRUE(overwrite)){
  rtracklayer::export(frag, filename_fh, format = "bam")
  message("bam file written successfully.")
} else {
  message("output file exists. \n")
  message(filename, " skipped. \n")
}



