
# Loading packages
library(httr)
library(jsonlite)
library(tidyverse)

# Initializing API Call
call <- "http://finaledb.research.cchmc.org/api/v1/seqrun?tissue=blood+plasma&assay=WGS&limit=10000"

# Getting details in API
get_movie_details <- GET(url = call)


# Converting content to text
get_movie_text <- content(get_movie_details,
                          "text", encoding = "UTF-8")

# Parsing data in JSON
get_movie_json <- fromJSON(get_movie_text,
                           flatten = TRUE)

# Converting into dataframe
get_movie_dataframe <- as.data.frame(get_movie_json)


files_list <- get_movie_dataframe$results.analysis.hg19

.helper <- function(x){
  return(paste0("https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/", x[["key"]]))
}

ans <- lapply(files_list, .helper)

ans_unlist <- unlist(ans)

.download_helper <- function(url){
  
  download.file(url = url, destfile = basename(url))
  
}

wget_helper <- function(url){
  
  system(paste("wget", url,  "--random-wait",  "-c", "--no-verbose"))
  
  
}

frag_tsv <- ans_unlist[stringr::str_detect(ans_unlist, "frag.tsv.bgz")]

job_log <- lapply(frag_tsv, wget_helper)




# meta

get_filename <- function(x){
  
  if("fragment" %in% x[["desc"]]) {
  ans <- dplyr::filter(x, desc == 'fragment') %>%
    dplyr::pull(key)
  return(ans)
    
  } else {
    
    return(NA_real_)
  }
}

raw_meta <- get_movie_dataframe



meta <- raw_meta %>% 
  as_tibble() 
  
  
hg19_tsv <- meta[["results.analysis.hg19"]]

hg19_tsv_filename <- lapply(hg19_tsv, get_filename) %>%
  unlist() %>%
  basename()
  
hg38_tsv <- meta[["results.analysis.hg38"]]

hg38_tsv_filename <- lapply(hg38_tsv, get_filename) %>%
  unlist() %>%
  basename()

meta_final <- meta %>% 
  dplyr::select(results.sample.id, 
                results.sample.name, 
                results.sample.age, 
                results.sample.gender, 
                results.sample.tissue, 
                results.sample.disease, 
                results.seqConfig.readlen, 
                results.seqConfig.seq_layout, 
                results.seqConfig.instrument,
                results.publication.author)


meta_final$frag_hg19 <- hg19_tsv_filename
meta_final$frag_hg38 <- hg38_tsv_filename



saveRDS(meta_final, file = "/Users/wang04/Documents/libs/github/ulyses/finaleDB/finaleDB_meta.rds")







