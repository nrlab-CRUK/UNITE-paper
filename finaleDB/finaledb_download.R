
# Loading packages
library(httr)
library(jsonlite)

setwd("/mnt/scratchc/nrlab/wang04/ulyses/data/external/finaledb")

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

frag_tsv <- ans_unlist[stringr::str_detect(ans_unlist, "frag.tsv.bgz")]


wget_helper <- function(url){
	  
	  system(paste("wget", url,  "--random-wait",  "-c", "--no-verbose"))
  
  
}


job_log <- lapply(frag_tsv, wget_helper)






