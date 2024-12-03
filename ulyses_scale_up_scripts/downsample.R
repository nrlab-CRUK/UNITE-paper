library(cfDNAPro)

args = commandArgs(trailingOnly=TRUE)

bamfile  <- as.character(args[1])

output_depth <- as.numeric(args[2])

output_dir <- as.character(args[3])

cfDNAPro::downsampleBam(bamfile = bamfile, 
                     genome = "hg19", 
                     output_depth = output_depth, 
		     output_type = c("bam"),
		     return_result = FALSE,
		     output_dir = output_dir,
                     flag = Rsamtools::scanBamFlag( isPaired = TRUE, 
                                                    isProperPair = NA, 
                                                    isUnmappedQuery = FALSE, 
                                                    hasUnmappedMate = FALSE, 
                                                    isMinusStrand = NA, 
                                                    isMateMinusStrand = NA,
                                                    isFirstMateRead = NA, 
                                                    isSecondMateRead = NA, 
                                                    isSecondaryAlignment = NA, 
                                                    isNotPassingQualityControls = NA,
                                                    isDuplicate = NA, 
                                                    isSupplementaryAlignment = NA),
                     what  = Rsamtools::scanBamWhat(),
                     tag = c("NM", "MD", "MC", "RG", "AS", "XS"))
