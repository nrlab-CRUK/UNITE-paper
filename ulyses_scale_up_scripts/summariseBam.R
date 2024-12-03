library(cfDNAPro)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

bamfile <- as.character(args[1])

output_dir <- as.character(args[2])

overwrite <- as.logical(args[3])



result <- cfDNAPro::summariseBam(
  bamfile = bamfile,
  total_count = TRUE,
  total_mapped_count = TRUE,
  chrM_count = TRUE,
  duplicate_count = TRUE,
  coverage_by_mapped = TRUE,
  gc_metrics = FALSE,
  loci_coverage_metrics = FALSE
)


filename <- file.path(output_dir, paste0(basename(bamfile), ".summary.csv"))


if ((!file.exists(filename)) || isTRUE(overwrite)) {
  write_csv(x = result, file = filename)
  message("csv file written successfully.")
} else {
  message("output file exists. \n")
  message(filename, " skipped. \n")
}
