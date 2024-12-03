
###Code for downloading data

```R

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

wget_helper <- function(url){
  
  system(paste("wget", url,  "--random-wait",  "-c", "--no-verbose"))
  
  
}

frag_tsv <- ans_unlist[stringr::str_detect(ans_unlist, "frag.tsv.bgz")]

job_log <- lapply(frag_tsv, wget_helper)



```

---

## unzip bgz files


```shell


for in in *.bgz
do
sbatch --time=0-3 --mem=16G -J "bgzip" --wrap="source ~/.bashrc; conda activate projectx; bgzip -d ${i}"
done

```

---

## convert tsv file to GRanges object in R

```R

chr_to_keep <- c(paste("chr", 1:22, sep = ""), "chrMT")

tmp <- read_csv("/Users/wang04/Documents/phd_projects/0_data_purge/finaleDB/EE85723.hg19.frag.tsv")

colnames(tmp) <- c("seqnames", "start", "end", "mapq", "strand")


tmp$seqnames <- paste0("chr", tmp[["seqnames"]])

tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(tmp, 
                                                  starts.in.df.are.0based = TRUE, 
                                                  keep.extra.columns = TRUE)


gr <- keepSeqlevels(tmp_gr, chr_to_keep , pruning.mode="coarse" )


seqinfo(gr) <- merge(seqinfo(gr), Seqinfo(genome = "hg19"))

gr <- keepSeqlevels(gr, chr_to_keep , pruning.mode="coarse" )


```

---

## meta data carpentry

```R
# convert data frame to tibble in r 
ans <- as_tibble()

```