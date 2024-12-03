


gatk <- function(gatk_path=gatk_path, method, args, maxheap="256G"){
  
  gatk_path <- "/scratcha/nrlab/wang04/urine/np_scripts/gatk-4.2.0.0/gatk"

  args <- cbind(paste("--java-options", paste0('"',"-Xmx", maxheap,'"'), method),
                args)
  retcode <- system2(gatk_path, args, stdout=FALSE)
  if(retcode != 0){
    stop( paste("GATK command [ java", paste(args, collapse=" "), "] failed."))
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
  overwrite_output = FALSE) {
                                
 if(!file.exists(input)) stop(paste("downsampleSamGATK cannot find", input))

 if (!overwrite_output) {
    if(file.exists(output)) stop("downsampleSamGATK failed to create", output, "file already exists.")
 } else {
    unlink(output)
 }

 if (as.numeric(depth) > 10) strategy <- "Chained"
  
  args <- paste("--INPUT", input,
                "--OUTPUT", output,
                "--STRATEGY", strategy,
                "--PROBABILITY", probability,
                "--ACCURACY", accuracy,
                "--CREATE_INDEX", create_index,
                "--METRICS_FILE", metrics_file,
                "--VALIDATION_STRINGENCY", stringency,
                "--VERBOSITY", verbosity,
                "--QUIET", quiet)

  gatk(method = "DownsampleSam", args = args)

}

# 'PostionBasedDownsampleSam' GATK function, this keeps the optical duplicates rates 
# steady, I feel this method is better than the 'DownsampleSam' function.
# see: https://gatk.broadinstitute.org/hc/en-us/articles/360036484732-PositionBasedDownsampleSam-Picard-
# and https://gatk.broadinstitute.org/hc/en-us/articles/360047717551-Downsampling-in-GATK

positionBasedDownsampleSamGATK <- function(input, 
  output, 
  fraction, 
  create_index = TRUE, 
  stringency = "LENIENT",
  quiet = FALSE,
  verbosity = "ERROR",
  overwrite_output = FALSE) {
                                
 if(!file.exists(input)) stop(paste("positionBasedDownsampleSamGATK cannot find", input))

 if (!overwrite_output) {
    if(file.exists(output)) stop("positonBasedDownsampleSamGATK failed to create", output, "file already exists.")
 } else {
    unlink(output)
 }

  args <- paste("--INPUT", input,
                "--OUTPUT", output,
                "--FRACTION", fraction,
                "--CREATE_INDEX", create_index,
                "--VALIDATION_STRINGENCY", stringency,
                "--VERBOSITY", verbosity,
                "--QUIET", quiet)

  gatk(method = "PositionBasedDownsampleSam", args = args)

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
  overwrite_output = FALSE) {
                                
 if(!file.exists(input)) stop(paste("markDuplicatesGATK cannot find", input))

 if (!overwrite_output) {
    if(file.exists(output)) stop("markDuplicatesGATK failed to create", output, "file already exists.")
 } else {
    unlink(output)
 }

  
  args <- paste("--INPUT", input,
                "--OUTPUT", output,
                "--CREATE_INDEX", create_index,
                "--METRICS_FILE", metrics_file,
                "--VERBOSITY", verbosity,
                "--VALIDATION_STRINGENCY", stringency,
                "--QUIET", quiet)

  gatk(method = "MarkDuplicates", args = args)

}

positionBasedDownsampleSamGATK_markDuplicatesGATK <- function(input, 
  output, 
  fraction, 
  mrkdup_input, 
  mrkdup_output) {
  
 positionBasedDownsampleSamGATK(input = input, output = output, fraction = fraction)
 markDuplicatesGATK(input = mrkdup_input, output = mrkdup_output)
 unlink(mrkdup_input)
 unlink(paste0(gsub(".bam", "", mrkdup_input), ".bai"))

}

downsampleSamGATK_markDuplicatesGATK <- function(input, 
  output, 
  probability, 
  depth,
  mrkdup_input, 
  mrkdup_output) {
  
 downsampleSamGATK(input = input, output = output, probability = probability, depth = depth)
 markDuplicatesGATK(input = mrkdup_input, output = mrkdup_output)
 unlink(mrkdup_input)
 unlink(paste0(gsub(".bam", "", mrkdup_input), ".bai"))

}


