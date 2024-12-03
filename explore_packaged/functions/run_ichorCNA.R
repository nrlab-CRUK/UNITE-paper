library(tidyverse, help, pos = 2, lib.loc = NULL)
library(rslurm, help, pos = 2, lib.loc = NULL)

# change to self control, raw_bams. (whatever)
run_ichorcna_default_ctrl <- function(input_bam, 
                                      wigfile, 
                                      ichorCNA_outdir,
                                      bam_basename,
                                      ichorCNA_script = "/home/nrlab/resources/ichorCNA/rubicon/ichorcna_rubicon_haichao.sh") {
  
  args <- paste(input_bam, wigfile, ichorCNA_outdir, bam_basename)
  system2(ichorCNA_script, args)
}

run_ichorcna_self_ctrl <- function(input_bam, 
                                   wigfile, 
                                   ichorCNA_outdir,
                                   bam_basename,
                                   ichorCNA_script = "/home/nrlab/resources/ichorCNA/rubicon/ichorcna_rubicon_haichao_self_control.sh") {
  
  args <- paste(input_bam, wigfile, ichorCNA_outdir, bam_basename)
  system2(ichorCNA_script, args)
}

