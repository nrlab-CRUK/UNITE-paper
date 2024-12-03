# arg parser --------------------------------------------------------------
library(argparse)

source("/home/nrlab/wang04/ulyses/models/cnn_xgboost_dataset_hyperparams.R")

# Create the argument parser
parser <- ArgumentParser(description = "Parameter Parser")

# Set the default values

authors_default <- authors_hyper 
cohorts_default <- cohorts_hyper
timepoints_default <- NULL


ulyses_dir_default <- "/home/nrlab/wang04/ulyses"
tensor_maker_code_default <- file.path(ulyses_dir_default, "explore_packaged", "functions/tensor_maker.r")
ichorcna_tf_cutoff_lower_default <- c(0)
ichorcna_tf_cutoff_upper_default <- c(1)
mae_file_default <- "/scratchc/nrlab/wang04/ulyses/data/MAE3.rds"
outdir_default <- "/Users/wang04/cnn_strat"
layers_default <- c("n_isize", "n_motif_smono1_C", "n_motif_smono1_T")
sample_type_default <- c("plasma")

# Add the arguments with default values
parser$add_argument("--ulyses_dir", dest = "ulyses_dir", default = ulyses_dir_default, help = "Directory path for ulyses project")
parser$add_argument("--tensor_maker_code", dest = "tensor_maker_code", default = tensor_maker_code_default, help = "File path for tensor_maker.r code")
parser$add_argument("--ichorcna_tf_cutoff_lower", dest = "ichorcna_tf_cutoff_lower", nargs = "+", type = "numeric", default = ichorcna_tf_cutoff_lower_default, help = "Lower cutoff values for ichorcna_tf_cutoff")
parser$add_argument("--ichorcna_tf_cutoff_upper", dest = "ichorcna_tf_cutoff_upper", nargs = "+", type = "numeric", default = ichorcna_tf_cutoff_upper_default, help = "Upper cutoff values for ichorcna_tf_cutoff")
parser$add_argument("--timepoints", dest = "timepoints", nargs = "*", default = timepoints_default, help = "List of timepoints")
parser$add_argument("--mae_file", dest = "mae_file", default = mae_file_default, help = "File path for MAE data")
parser$add_argument("--outdir", dest = "outdir", default = outdir_default, help = "Output directory path")
parser$add_argument("--layers", dest = "layers", nargs = "+", default = layers_default, help = "Names of layers")
parser$add_argument("--authors", dest = "authors", nargs = "+", default = authors_default, help = "Names of authors")
parser$add_argument("--cohorts", dest = "cohorts", nargs = "+", default = cohorts_default, help = "Names of cohorts")
parser$add_argument("--sample_type", dest = "sample_type", default = sample_type_default, help = "String vector, Sample type")


# Parse the arguments
args <- parser$parse_args()

# Access the parsed arguments
ulyses_dir <- args$ulyses_dir
tensor_maker_code <- args$tensor_maker_code
ichorcna_tf_cutoff_lower <- args$ichorcna_tf_cutoff_lower
ichorcna_tf_cutoff_upper <- args$ichorcna_tf_cutoff_upper
timepoints <- args$timepoints
mae_file <- args$mae_file
outdir <- args$outdir
layers <- args$layers
authors <- args$authors
cohorts <- args$cohorts
sample_type <- args$sample_type

source(tensor_maker_code)


# create a tibble to store parameters
para_df <- tibble(
                  ichorcna_tf_cutoff_lower = ichorcna_tf_cutoff_lower, 
                  ichorcna_tf_cutoff_upper = ichorcna_tf_cutoff_upper)

# run tensor_maker using para_df as input
para_df %>% 
  furrr::future_pmap(.f = tensor_maker,
  .progress = TRUE,
  layers = layers,
  authors = authors,
  cohorts = cohorts,
  sample_type = sample_type,
  timepoints = timepoints,
  mae_file = mae_file,
  outdir = outdir
  )


message("Done!")



