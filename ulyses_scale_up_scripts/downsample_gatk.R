args <- commandArgs(trailingOnly = TRUE)

bam_stats <- args[1]
target_depth <- as.numeric(args[2])
out_dir <- args[3]
gatk_path <- args[4]
random_or_pos <- args[5]

if (is.null(gatk_path) || is.na(gatk_path)) {
  gatk_path <- "/home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"
}

if (is.null(random_or_pos) || is.na(random_or_pos)) {
  random_or_pos <- "random"
}

setwd(dirname(bam_stats))
library(dplyr)
library(readr)

bam_stats <- read_csv(bam_stats)

if (bam_stats[["coverage_by_mapped_reads"]] < target_depth) {
  stop("no enought reads for ", target_depth)
}
# if random_or_pos is random

if (random_or_pos == "random") {
  params <- bam_stats %>%
    dplyr::filter(coverage_by_mapped_reads >= target_depth) %>%
    dplyr::mutate(input = file) %>%
    dplyr::mutate(depth = coverage_by_mapped_reads) %>%
    dplyr::mutate(probability = target_depth / depth) %>%
    dplyr::mutate(output = file.path(!!out_dir, paste0(basename(input), ".", target_depth, "x", ".bam"))) %>%
    dplyr::mutate(mrkdup_input = output) %>%
    dplyr::mutate(mrkdup_output = file.path(out_dir, paste0(basename(input), ".", target_depth, "x", ".mrkdup.bam"))) %>%
    dplyr::select(input, output, probability, depth, mrkdup_input, mrkdup_output)
} else if (random_or_pos == "pos") {
  params <- bam_stats %>%
    dplyr::filter(coverage_by_mapped_reads >= target_depth) %>%
    dplyr::mutate(input = file) %>%
    dplyr::mutate(depth = coverage_by_mapped_reads) %>%
    dplyr::mutate(probability = target_depth / depth) %>%
    dplyr::mutate(output = file.path(!!out_dir, paste0(basename(input), ".", target_depth, "x", ".posbased.bam"))) %>%
    dplyr::mutate(mrkdup_input = output) %>%
    dplyr::mutate(mrkdup_output = file.path(out_dir, paste0(basename(input), ".", target_depth, "x", ".posbased.mrkdup.bam"))) %>%
    dplyr::select(input, output, probability, depth, mrkdup_input, mrkdup_output)
} else {
  stop("random_or_pos must be either random or pos")
}





input <- params[["input"]]
output <- params[["output"]]
probability <- params[["probability"]]
depth <- params[["depth"]]
mrkdup_input <- params[["mrkdup_input"]]
mrkdup_output <- params[["mrkdup_output"]]






gatk <- function(gatk_path,
                 method,
                 args,
                 maxheap = "300G") {
  args <- cbind(
    paste("--java-options", paste0('"', "-Xmx", maxheap, '"'), method),
    args
  )
  retcode <- system2(gatk_path, args, stdout = "")
  if (retcode != 0) {
    stop(paste("GATK command [ java", paste(args, collapse = " "), "] failed."))
  }
  retcode
}

downsampleSamGATK <- function(input,
                              output,
                              probability,
                              depth,
                              strategy = "HighAccuracy",
                              accuracy = 0.00001,
                              create_index = TRUE,
                              metrics_file = paste0(output, ".metrics"),
                              quiet = FALSE,
                              stringency = "LENIENT",
                              verbosity = "ERROR",
                              overwrite_output = FALSE, ...) {
  if (!file.exists(input)) stop(paste("downsampleSamGATK cannot find", input))

  if (!overwrite_output) {
    if (file.exists(output)) stop("downsampleSamGATK failed to create", output, "file already exists.")
  } else {
    unlink(output)
  }

  if (as.numeric(depth) > 10) strategy <- "Chained"

  args <- paste(
    "--INPUT", input,
    "--OUTPUT", output,
    "--STRATEGY", strategy,
    "--PROBABILITY", probability,
    "--ACCURACY", accuracy,
    "--CREATE_INDEX", create_index,
    "--METRICS_FILE", metrics_file,
    "--VALIDATION_STRINGENCY", stringency,
    "--VERBOSITY", verbosity,
    "--QUIET", quiet
  )

  gatk(method = "DownsampleSam", args = args, ...)
}


positionBasedDownsampleSamGATK <- function(input,
                                           output,
                                           fraction,
                                           create_index = TRUE,
                                           stringency = "LENIENT",
                                           quiet = FALSE,
                                           verbosity = "ERROR",
                                           overwrite_output = FALSE, ...) {
  if (!file.exists(input)) stop(paste("positionBasedDownsampleSamGATK cannot find", input))

  if (!overwrite_output) {
    if (file.exists(output)) stop("positonBasedDownsampleSamGATK failed to create", output, "file already exists.")
  } else {
    unlink(output)
  }

  args <- paste(
    "--INPUT", input,
    "--OUTPUT", output,
    "--FRACTION", fraction,
    "--CREATE_INDEX", create_index,
    "--VALIDATION_STRINGENCY", stringency,
    "--VERBOSITY", verbosity,
    "--QUIET", quiet
  )

  gatk(method = "PositionBasedDownsampleSam", args = args, ...)
}

# markDuplicatesGATK function
# markDuplicatesGATK function
markDuplicatesGATK <- function(input,
                               output,
                               create_index = TRUE,
                               metrics_file = paste0(output, ".metrics"),
                               quiet = FALSE,
                               stringency = "LENIENT",
                               verbosity = "ERROR",
                               overwrite_output = FALSE, ...) {
  if (!file.exists(input)) stop(paste("markDuplicatesGATK cannot find", input))

  if (!overwrite_output) {
    if (file.exists(output)) stop("markDuplicatesGATK failed to create", output, "file already exists.")
  } else {
    unlink(output)
  }


  args <- paste(
    "--INPUT", input,
    "--OUTPUT", output,
    "--CREATE_INDEX", create_index,
    "--METRICS_FILE", metrics_file,
    "--VERBOSITY", verbosity,
    "--VALIDATION_STRINGENCY", stringency,
    "--QUIET", quiet
  )

  gatk(method = "MarkDuplicates", args = args, ...)
}

positionBasedDownsampleSamGATK_markDuplicatesGATK <- function(input,
                                                              output,
                                                              fraction,
                                                              mrkdup_input,
                                                              mrkdup_output, ...) {
  positionBasedDownsampleSamGATK(input = input, output = output, fraction = fraction, ...)
  markDuplicatesGATK(input = mrkdup_input, output = mrkdup_output, ...)
  unlink(mrkdup_input)
  unlink(paste0(gsub(".bam", "", mrkdup_input), ".bai"))
}

downsampleSamGATK_markDuplicatesGATK <- function(input,
                                                 output,
                                                 probability,
                                                 depth,
                                                 mrkdup_input,
                                                 mrkdup_output, ...) {
  downsampleSamGATK(input = input, output = output, probability = probability, depth = depth, ...)
  markDuplicatesGATK(input = mrkdup_input, output = mrkdup_output, ...)
  unlink(paste0(gsub(".bam$", "", mrkdup_input), ".bai"))
  unlink(mrkdup_input)
  unlink(paste0(mrkdup_input, ".metrics"))
}


if (random_or_pos == "random") {
  message("downsampleSamGATK_markDuplicatesGATK...")
  downsampleSamGATK_markDuplicatesGATK(
    input = input,
    output = output,
    probability = probability,
    depth = depth,
    mrkdup_input = mrkdup_input,
    mrkdup_output = mrkdup_output,
    gatk_path = gatk_path
  )
} else if (random_or_pos == "pos") {
  message("positionBasedDownsampleSamGATK_markDuplicatesGATK...")
  positionBasedDownsampleSamGATK_markDuplicatesGATK(
    input = input,
    output = output,
    fraction = probability,
    mrkdup_input = mrkdup_input,
    mrkdup_output = mrkdup_output,
    gatk_path = gatk_path
  )
}
