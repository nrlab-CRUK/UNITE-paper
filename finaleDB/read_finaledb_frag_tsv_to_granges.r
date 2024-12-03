
library(GenomicAlignments)
library(GenomeInfoDb)

args <- commandArgs(trailingOnly = TRUE)

#-----------------------------------------------------------------------------
#params
input <- args[1]
genome_label <- args[2] 

if(is.null(genome_label) || is.na(genome_label)){

        if(stringr::str_detect(input, "hg19")){
                genome_label <- "hg19"
        } else if(stringr::str_detect(input, "hg38")){
                genome_label <- "hg38"
        }

}

output <- paste(input, ".GRanges.rds", sep = "")
#-----------------------------------------------------------------------------

tmp <- readr::read_tsv(input, col_names = FALSE, show_col_types = FALSE)

colnames(tmp) <- c("seqnames", "start", "end", "mapq", "strand")

if(!stringr::str_detect(tmp$seqnames[[1]], "^chr")){
        tmp$seqnames <- paste0("chr", tmp[["seqnames"]])
}


tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(tmp, 
                                                  starts.in.df.are.0based = TRUE, 
                                                  keep.extra.columns = TRUE)


seqinfo(tmp_gr) <- S4Vectors::merge(seqinfo(tmp_gr), Seqinfo(genome = genome_label))

gr <- keepStandardChromosomes(tmp_gr,pruning.mode = "coarse", species = "Homo sapiens")


saveRDS(gr, file = output)

message("saved ", output, "\n" )

