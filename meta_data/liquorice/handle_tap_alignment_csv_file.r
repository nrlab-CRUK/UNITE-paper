

dt <- read_csv("/Users/wang04/Documents/phd_projects/0_data_purge/external/liquorice/liquorice_tap_alignment.csv")

write_csv(x = dt, file = "/Users/wang04/Documents/phd_projects/0_data_purge/external/liquorice/liquorice_tap_alignment_quoted.csv",na = "", quote = "all")
