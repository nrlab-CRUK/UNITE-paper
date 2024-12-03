
library(tidyverse)

dt <- read_csv("/Users/wang04/Documents/phd_projects/0_data_purge/external/delfi/delfi_tap_alignment.csv")

dt <- dt %>% dplyr::mutate(SequencingDate = "2023-01-01")

write_csv(x = dt, file = "/Users/wang04/Documents/phd_projects/0_data_purge/external/delfi/delfi_tap_alignment_quoted.csv",na = "", quote = "all")
