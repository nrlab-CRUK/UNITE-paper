
source("/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/frag_functions.r")


    # FUN to get the N mismatch on positive and negative strands
    
    get_pos_nm <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            mcols(GenomicAlignments::first(galp_object))$NM,
            mcols(GenomicAlignments::second(galp_object))$NM
        )
    }
        
    get_neg_nm <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            mcols(GenomicAlignments::first(galp_object))$NM,
            mcols(GenomicAlignments::second(galp_object))$NM
        )
    }
    
    # FUN to get the map quality of positive and negative strands
    
    get_pos_mapq <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            mcols(GenomicAlignments::first(galp_object))$mapq,
            mcols(GenomicAlignments::second(galp_object))$mapq
        )
    }
        
    get_neg_mapq <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            mcols(GenomicAlignments::first(galp_object))$mapq,
            mcols(GenomicAlignments::second(galp_object))$mapq
        )
    }
    

    # FUN to get the width of positive and negative strands
    
    get_pos_width <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            width(GenomicAlignments::first(galp_object)),
            width(GenomicAlignments::second(galp_object))
        )
    }
    
    
    get_neg_width <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            width(GenomicAlignments::first(galp_object)),
            width(GenomicAlignments::second(galp_object))
        )
    }
      
    # FUN to get the CIGAR of positive and negative strands
    
    get_pos_cigar <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            GenomicAlignments::cigar(GenomicAlignments::first(galp_object)),
            GenomicAlignments::cigar(GenomicAlignments::second(galp_object))
        )
    }
    
    
    get_neg_cigar <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            GenomicAlignments::cigar(GenomicAlignments::first(galp_object)),
            GenomicAlignments::cigar(GenomicAlignments::second(galp_object))
        )
    }
    


bam_to_galp_proper_pair <- function(bamfile,
                        chromosome_to_keep = paste0("chr", 1:22),
                        strand_mode = 1,
                        ...) {
    # Check parameters
    
    #if (!file.exists(bamfile)) {
    #    m <- paste0(bamfile, " doesn't exist, please check your input.")
    #    stop(m)
    #}
    
    
    # Read bam into galp
    
    flag <- Rsamtools::scanBamFlag(
        isPaired = TRUE,
 		 isProperPair = TRUE,
        isDuplicate = FALSE,
        isSecondaryAlignment = FALSE,
        isUnmappedQuery = FALSE,
        isSupplementaryAlignment = FALSE
    )
    
    what <- c("cigar", "mapq", "isize")
    
    # MD: string encoding mismatched and deleted reference bases
    # NM: number of differences (mismatch + I + D) between the seq and ref
    # XM: Number of mismatches in the alignment (specific to bwa)
    # XO: Number of gap opens (specific to bwa)
    # XG: Number of gap extensions (specific to bwa)
    
    tag <- c("MD", "NM", "XM", "XO", "XG")
    
    
    param <- Rsamtools::ScanBamParam(flag = flag,
                                     what = what,
                                     tag = tag)
    
    # read bam file as galp object
    galp <- GenomicAlignments::readGAlignmentPairs(
        bamfile,
        index = bamfile,
        param = param,
        use.names = TRUE,
        strandMode = strand_mode
    )
    
    # strandMode should be one for downstream operations
    stopifnot(GenomicAlignments::strandMode(galp) == 1)
    
    # remove read pairs without correct seqnames and strand information
    
    message("Curating seqnames and strand information...")
    # only keep chr1:22
    galp <-
        keepSeqlevels(galp, chromosome_to_keep, pruning.mode = "coarse")
    
    # remove read pairs without seqname information
    galp2 <- galp[!is.na(GenomicAlignments::seqnames(galp))]
    
    # remove read pairs without strand information
    galp3 <- galp2[strand(galp2) != '*']
    
    return(galp3)
    
}







bam_to_frag <- function(
x,
genome = "hg19",
filter_by_proper_pair = TRUE,
curate_start_end = TRUE,
isize_from = 1,
isize_to = 1000) {
	library(GenomicAlignments)

	if (genome == "hg19") {
	library(BSgenome.Hsapiens.UCSC.hg19)
	genome <- BSgenome.Hsapiens.UCSC.hg19
	} else if (genome == "hg38") {
	library(BSgenome.Hsapiens.UCSC.hg38)
	genome <- BSgenome.Hsapiens.UCSC.hg38
	}

	if (filter_by_proper_pair) {
	galp <- bam_to_galp_proper_pair(x) 
	} else {
	galp <- bam_to_galp(x) 
	}

	galp <- remove_outward_facing_readpairs(galp)


	if (curate_start_end) {

		frag <- curate_start_and_end(galp)

	} else {

		frag <- granges(galp)

       frag$pos_nm <- get_pos_nm(galp)
       frag$neg_nm <- get_neg_nm(galp)
       frag$pos_mapq <- get_pos_mapq(galp)
       frag$neg_mapq <- get_neg_mapq(galp)
       frag$pos_width <- get_pos_width(galp)
       frag$neg_width <- get_neg_width(galp)
       frag$pos_cigar <- get_pos_cigar(galp)
       frag$neg_cigar <- get_neg_cigar(galp)
		frag$frag_len  <- BiocGenerics::width(frag)
		frag$frag_len_seen <- 1L

	}

	frag <- frag %>%
    		`seqlengths<-`(seqlengths(genome)[levels(seqnames(.))]) %>%
    		`genome<-`(seqinfo(genome)@genome[1])


	frag <- remove_out_of_bound_reads(frag) %>%
			plyranges::filter(frag_len >= isize_from & frag_len <= isize_to) %>%
			plyranges::filter(pos_mapq >= 0 & neg_mapq >= 0) %>%
			plyranges::filter(!stringr::str_detect(pos_cigar, pattern = "[ID]")) %>%
			plyranges::filter(!stringr::str_detect(neg_cigar, pattern = "[ID]")) 

  
  return(frag)

}


read_isize <- function(x,
	genome = "hg19",
	filter_by_proper_pair = TRUE,
	curate_start_end = TRUE,
	isize_from = 1,
	isize_to = 1000){

	frag <- bam_to_frag(x,
		genome = genome,
		filter_by_proper_pair = filter_by_proper_pair ,
		curate_start_end = curate_start_end,
		isize_from = isize_from,
		isize_to = isize_to)
  
  isize_tibble <- tibble("insert_size" = frag$frag_len , "count" = frag$frag_len_seen )
  
  result <- isize_tibble %>%
    dplyr::group_by(.data$insert_size) %>%
    dplyr::summarise("All_Reads.fr_count" = sum(count))
  
  return(result)

}







isize_plot <- function(x){

	df <- read_isize(x)
	message("read bams done!")

	p <- ggplot(df) + 
		geom_line(aes(insert_size, All_Reads.fr_count)) + 
		geom_vline(xintercept = c(167, 333), linetype = "dashed", color = "darkgrey") +
  		geom_ribbon(aes(x = insert_size, ymax = All_Reads.fr_count), ymin = 0, alpha=0.3) +
  		scale_fill_manual(name='', values = "grey") +
  		xlab("Fragment Length (bp)") +
  		ylab("Count") +
		coord_cartesian(xlim = c(30, 500))+ 
  		scale_x_continuous(limits = c(20, 500), 
                     breaks = c(30, 100, 167, 200, 333, 400, 500), 
                     labels = c("30",  "100", "167", "200", "333", "400", "500")) + 
  		theme_classic()

	message("saving pdf files:")
	message(paste0(x, ".refined_isize_dist.pdf"))
	ggsave(filename = paste0(x, ".refined_isize_dist.pdf"), plot = p, width = 8, height = 5)

	saveRDS(object = df, paste0(x, ".isize_df.rds"))
	message("saved df rds file, job done")
	return(df)
}
