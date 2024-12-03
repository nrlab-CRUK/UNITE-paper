library(readr)
library(tibble)


ref_out_dir = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/resources"

ref <- tibble(output_col = c("n_isize", 
                             "n_pos_nm", 
                             "n_neg_nm",
                             "n_motif_s1_A",
                             "n_motif_s1_T",
                             "n_motif_s1_C",
                             "n_motif_s1_G",
                             "bin_mean_GC",
                             "bin_mean_mappability"),
              
              source_col = c("frag_len_seen", 
                             "pos_nm",
                             "neg_nm",
                             "motif_s1",
                             "motif_s1",
                             "motif_s1",
                             "motif_s1",
                             "gc",
                             "mappability"),
              
              method  = c("sum",
                          "sum",
                          "sum",
                          "motif_filter_nrow",
                          "motif_filter_nrow",
                          "motif_filter_nrow",
                          "motif_filter_nrow",
                          "mean",
                          "mean"
                          ),
              
              null_filling_value = c(0, 
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     0,
                                     NA_real_,
                                     NA_real_
                                     ),
              type = c("length",
                       "mismatch",
                       "mismatch",
                       "motif",
                       "motif",
                       "motif",
                       "motif",
                       "bin_status",
                       "bin_status"),
              
              type_for_norm = c("length",
                       "mismatch",
                       "mismatch",
                       "motif_s1",
                       "motif_s1",
                       "motif_s1",
                       "motif_s1",
                       "bin_status_GC",
                       "bin_status_mappability")
              
              )

ref_filename <- file.path(ref_out_dir, "overlap_ref.csv")
readr::write_csv(ref, ref_filename)