hostname <- 'inst-webappdev.cri.camres.org'
username <- 'nrlab_reader'
password <- 'maryhadalittleiguana'
database <- 'wordpress'

library(tidyverse)

con <- DBI::dbConnect(RMariaDB::MariaDB(),
                      dbname = database,
                      user = username,
                      password = password,
                      host = hostname)

DBI::dbListTables(con)


wd <- "/Users/wang04/Documents/libs/github/ulyses/meta_data/DBS"
#wd <- "/home/nrlab/wang04/ulyses/meta_data/DBS"
outfile <- file.path(wd, "dbs_meta.csv")


slxid <- c("SLX-20706", 
           "SLX-22151", 
           "SLX-21409", 
           "SLX-22439", 
           "SLX-22440", 
           "SLX-19894", 
           "SLX-21959", 
           "SLX-21619")

slxid2 <- c("SLX-20706", 
           "SLX-22151", 
           "SLX-21409", 
           "SLX-22439", 
           "SLX-22440", 
           "SLX-19894", 
           "SLX-21959", 
           "SLX-21619",
           "SLX-22916",
           "SLX-22917")

expid <- c("EXP3161_AA",
           "EXP3221_AA",
           "EXP3179_AA",
           "EXP3235_AA",
           "EXP3235_AA",
           "EXP3123_AA",
           "EXP3200",
           "EXP3184")




get_meta <- function(expid, slxid){
  
  
  cases_data <- tbl(con, "rosenfeld_DNASeq_layout") %>% 
    filter(expid == expid ) %>% 
    filter(SLX_ID == slxid) %>% 
    left_join(tbl(con, "rosenfeld_Sample"), by = c("sampleName" = "sampleName")) %>% 
    select(sampleName, SLX_ID, barcode, clinicalCaseName) %>% 
    left_join(tbl(con, "rosenfeld_Cases"), by = c("clinicalCaseName" = "caseName")) %>% 
    collect()
  
  
  return(cases_data)
  
  
  
}


get_meta2 <- function(slxid){
  
  
  cases_data <- tbl(con, "rosenfeld_DNASeq_layout") %>% 
    filter(SLX_ID == slxid) %>% 
    left_join(tbl(con, "rosenfeld_Sample"), by = c("sampleName" = "sampleName")) %>% 
    select(sampleName, SLX_ID, barcode, clinicalCaseName) %>% 
    left_join(tbl(con, "rosenfeld_Cases"), by = c("clinicalCaseName" = "caseName")) %>% 
    collect()
  
  
  return(cases_data)
  
  
  
}

#ans <- map2(expid, slxid, get_meta) %>%
  list_rbind()

ans <- map(slxid, get_meta2) %>%
  list_rbind()


ans <- ans %>%
  # only keep SLX-19894 (only D711tp_D502tp, D711tp_D504tp)
  mutate(keep = case_when(
    SLX_ID == "SLX-19894" & barcode %in% c("D711tp-D502tp", "D711tp-D504tp") ~ 1,
    SLX_ID %in% c("SLX-20706", 
                  "SLX-22151", 
                  "SLX-21409", 
                  "SLX-22439", 
                  "SLX-22440", 
                  "SLX-21959", 
                  "SLX-21619",
                  "SLX-22916",
                  "SLX-22917") ~ 1,
    TRUE ~ 0
  )) %>%
  filter(keep == 1) %>%
  add_row(SLX_ID = c("SLX-22151", "SLX-22439", "SLX-22440"))

write_csv(x = ans, file = outfile)


############





