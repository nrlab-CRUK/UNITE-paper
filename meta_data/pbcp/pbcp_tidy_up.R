library(tidyverse)
library(readxl)
wd <- "/home/nrlab/wang04/ulyses/meta_data/pbcp"

archive_path <- "/archive/Groups/NRLab/fs07/Manuscripts/202305_PBCP_Santonja/Data/sWGS"

meta_xlsx <- "/home/nrlab/wang04/ulyses/meta_data/pbcp/PBCP_libraries_sWGS.xlsx"
id_xlsx <- "/home/nrlab/wang04/ulyses/meta_data/pbcp/clinical_data/Sex ethnicity data of patients_with publication ID.xlsx"

id <- readxl::read_xlsx(id_xlsx)

meta_raw <- readxl::read_xlsx(meta_xlsx )

meta_tidy <- meta_raw %>%
	mutate(bam_id = paste(`sWGS run in SLX`, `Index`, sep = ".")) %>%
	select(`Input_sample`, `Case/control`, bam_id) 

meta_tidy2 <- meta_tidy %>%
left_join(id, by = c("Input_sample" = "Nrlab Database ID NOT TO BE PUBLISHED"))

# write meta_tidy2 as csv
write_csv(meta_tidy2, file.path(wd, "pbcp_meta_tidy.csv"))
