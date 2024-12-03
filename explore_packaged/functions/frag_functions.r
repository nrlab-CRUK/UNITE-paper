
bam_to_galp <- function(bamfile,
                        chromosome_to_keep = paste0("chr", 1:22),
                        strand_mode = 1,
                        #which_genomic_region = which_genomic_region,
                        ...) {
    # Check parameters
    
    #if (!file.exists(bamfile)) {
    #    m <- paste0(bamfile, " doesn't exist, please check your input.")
    #    stop(m)
    #}
    
    
    # Read bam into galp
    
    flag <- Rsamtools::scanBamFlag(
        isPaired = TRUE,
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
                                     #which = which_genomic_region,
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




remove_outward_facing_readpairs <- function(galp) {
    # Get FR read pairs(i.e. first read aligned to positive strand)
    
    FR_id <-
        BiocGenerics::which(strand(GenomicAlignments::first(galp)) == "+")
    FR_galp <- galp[FR_id]
    
    # Get RF read pairs(i.e. second read aligned to positive strand)
    RF_id <-
        BiocGenerics::which(strand(GenomicAlignments::first(galp)) == "-")
    RF_galp <- galp[RF_id]
    
    # Get FR inward facing pairs
    FR_galp_inward <-
        FR_galp[GenomicAlignments::start(GenomicAlignments::first(FR_galp)) <
                    GenomicAlignments::end(GenomicAlignments::second(FR_galp))]
    
    # Get RF inward facing pairs
    RF_galp_inward <-
        RF_galp[GenomicAlignments::start(GenomicAlignments::second(RF_galp)) <
                    GenomicAlignments::end(GenomicAlignments::first(RF_galp))]
    
    message("Removing outward facing frags ...")
    # Combine RF and FR inward read pairs
    galp_corrected <- c(FR_galp_inward, RF_galp_inward)
    
    return(galp_corrected)
    
}




curate_start_and_end <- function(galp)  {
  
    message("Corrected start and end coordinates of fragments ...")
    # FUN to get the start and end of positive and negative strands
  
    get_start <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            start(GenomicAlignments::first(galp_object)),
            start(GenomicAlignments::second(galp_object))
        )
    }
    
    get_end <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            end(GenomicAlignments::first(galp_object)),
            end(GenomicAlignments::second(galp_object))
        )
    }

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
    

    # Add the fragment information back to fragmentwise object
    
    obj <-  GRanges(
        seqnames = seqnames(GenomicAlignments::first(galp)),
        ranges = IRanges(start = get_start(galp),
                         end = get_end(galp)),
        strand = strand(galp),
        pos_nm = get_pos_nm(galp),
        neg_nm = get_neg_nm(galp),
        pos_mapq = get_pos_mapq(galp),
        neg_mapq = get_neg_mapq(galp),
        pos_width = get_pos_width(galp),
        neg_width = get_neg_width(galp),
        pos_cigar = get_pos_cigar(galp),
        neg_cigar = get_neg_cigar(galp)
    )
    
    # Add seqnames
    
     names(obj) <- names(galp)
    
    
    # Add frag_len meta col
  
    obj <- obj %>%
      plyranges::mutate(frag_len = BiocGenerics::width(obj), frag_len_seen = 1) 

    return(obj)
}



get_out_of_bound_index <- function (x) {
    if (length(x) == 0L)
        return(integer(0))
    
    x_seqnames_id <- as.integer(seqnames(x))
    x_seqlengths <- unname(seqlengths(x))
    seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
    seqlength_is_na <- is.na(x_seqlengths)
    seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
    
    which(seqlevel_has_bounds[x_seqnames_id] &
              (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}

remove_out_of_bound_reads <- function(granges_object){
    idx <-get_out_of_bound_index(granges_object)
    message("Removing out-of-bound reads ...")
    if(length(idx) != 0L)  granges_object <- granges_object[-idx]
    return(granges_object)
}

# Motif annotation

get_motif <- function(obj,
                	genome,
                	motif_type, 
                	motif_length){

      obj <- unstrand(obj)

      if(motif_type == "s") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::start(obj),
                               end = GenomicAlignments::start(obj) + motif_length - 1))
      }

      # add a new motif type called 'smono' for single base motif annotation
        if(motif_type == "smono") {
        coord <- GRanges(seqnames = seqnames(obj),
                ranges = IRanges(start = GenomicAlignments::start(obj) + motif_length - 1,
                                 end = GenomicAlignments::start(obj) + motif_length - 1))
        }


      if(motif_type == "e") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::end(obj) - motif_length + 1 ,
                               end = GenomicAlignments::end(obj)))
      }

      if(motif_type == "emono") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::end(obj) - motif_length + 1 ,
                               end = GenomicAlignments::end(obj) - motif_length + 1))
      }


      if(motif_type == "u") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::start(obj) - motif_length,
                               end = GenomicAlignments::start(obj) - 1))
      }

      if(motif_type == "umono") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::start(obj) - motif_length,
                               end = GenomicAlignments::start(obj) - motif_length))
      }

      if(motif_type == "d") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::end(obj) + 1,
                               end = GenomicAlignments::end(obj) + motif_length))
      }

      if(motif_type == "dmono") {
      coord <- GRanges(seqnames = seqnames(obj),
              ranges = IRanges(start = GenomicAlignments::end(obj) + motif_length,
                               end = GenomicAlignments::end(obj) + motif_length))
      }
      
      seqlengths(coord) <- seqlengths(genome)[levels(seqnames(coord))]
      genome(coord) <- seqinfo(genome)@genome %>% unique()
        
      id <- get_out_of_bound_index(coord)
      
      if(length(id) != 0)  {
          coord <- plyranges::mutate(coord, index = 1:length(coord)) 
          coord_tidy <- coord[-id]
          coord_tidy <- plyranges::mutate(coord_tidy, anno  = getSeq(genome, coord_tidy, as.character = TRUE))
          df <- tibble::tibble(index = coord_tidy$index, anno = coord_tidy$anno)
          result <- dplyr::bind_rows(df, tibble(index = id)) %>% 
            dplyr::arrange(index) %>%
            dplyr::pull(anno)
        
      } else {
        result <- BSgenome::getSeq(genome, coord, as.character = TRUE)
      }
        
      return(result)
      
}

annotate_motif <- function(frag, 
                          genome,
                          motif_type_vec = factor(c("u", "s", "e", "d"), levels = c("u", "s", "e", "d")),
                          motif_length_vec = 1:3, ...){

                          
    motif_type_vec <- sort(motif_type_vec)


    for (motif_type in motif_type_vec ) {
      for(motif_length in motif_length_vec) {
        
        col_name <- paste0("motif", "_", motif_type, motif_length)

        message(paste0("Handling ", col_name, " ..." ))

        frag_new <- frag %>%
                    plyranges::mutate({{col_name}} := get_motif(obj = frag, 
                                                    genome = genome, 
                                                    motif_length = motif_length, 
                                                    motif_type = motif_type))
        frag <- frag_new

      }
    }

    return(frag)

}


# Filter frag

report_nrow <- function(obj, filter_step){

  n_row <- length(obj)

  msg <- paste0(filter_step, " ", n_row, " fragments left.")

  print(msg)

}





###############################################################################
# Simplified functions for compatability purpose
###############################################################################

curate_start_and_end_simplified <- function(galp)  {
  
    # FUN to get the start and end of positive and negative strands
  
    get_start <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "+",
            start(GenomicAlignments::first(galp_object)),
            start(GenomicAlignments::second(galp_object))
        )
    }
    
    get_end <- function(galp_object) {
        ifelse(
            strand(GenomicAlignments::first(galp_object)) == "-",
            end(GenomicAlignments::first(galp_object)),
            end(GenomicAlignments::second(galp_object))
        )
    }

    

    # Add the fragment information back to fragmentwise object
    
    obj <-  GRanges(
        seqnames = seqnames(GenomicAlignments::first(galp)),
        ranges = IRanges(start = get_start(galp),
                         end = get_end(galp)),
        strand = strand(galp)
    )
    
    message("Corrected start and end coordinates of fragments ...")

    # Add seqnames
     names(obj) <- names(galp)

    return(obj)
}




###############################################################################
# export RDS files new functions
###############################################################################



export_summary_table_long <- function(x, 
                                      layers){
  
  
  tibble_wide <-  dplyr::bind_rows(x, .id = "tile_size") %>%
    dplyr::mutate(unique_group_id = row_number()) %>%
    dplyr::mutate(chr = stringr::str_extract(id, pattern = "\\b\\d{1,2}") %>% as.integer()) %>%
    dplyr::mutate(tile_size_label = stringr::str_extract(tile_size, pattern = "\\d+\\b") %>% as.integer()) %>%
    tidyr::separate(col = id, into = c("arm", "within_arm_index"), sep = "_", remove = FALSE, convert = TRUE) %>%
    dplyr::relocate(chr, .after = arm) 
  
  tibble_wide$arm <- factor(tibble_wide$arm, levels = tibble_wide$arm %>% unique()) 
  tile_size_levels <- paste0("n_100kb_bins_", tibble_wide %>% dplyr::pull(tile_size_label) %>% sort() %>% unique() )  
  tibble_wide$tile_size <-  factor(tibble_wide$tile_size, levels = tile_size_levels)
  
  # arrange columns in specific order
  tibble_wide <- tibble_wide %>%
    dplyr::arrange(tile_size, id, within_arm_index, frag_len) %>%
    # careful
    dplyr::group_by(tile_size, frag_len) %>% 
    dplyr::mutate(id_rename = row_number() ) %>% 
    dplyr::relocate(id_rename, .after = id) %>%
    dplyr::relocate(tile_size_label, .after = tile_size) %>%
    dplyr::ungroup()
  
  tibble_long <- tibble_wide %>%
    tidyr::pivot_longer( 
      cols = all_of(layers), 
      names_to = "channel", 
      values_to = "pixel") %>%
    dplyr::mutate(unique_group_id = row_number())
  
  tibble_long$channel <- factor(tibble_long$channel , levels = layers)
  
  tibble_long <- tibble_long %>%
    # the order of cols in the arrange correspond to the 
    # dimension of ndarray in python analysis
    dplyr::arrange(tile_size, channel, id, frag_len)

  return(tibble_long)
  
}



export_summary_table_wide <- function(x){
  
  tibble_wide <-  dplyr::bind_rows(x, .id = "tile_size") %>%
    dplyr::mutate(unique_group_id = row_number()) %>%
    dplyr::mutate(chr = stringr::str_extract(id, pattern = "\\b\\d{1,2}") %>% as.integer()) %>%
    dplyr::mutate(tile_size_label = stringr::str_extract(tile_size, pattern = "\\d+\\b") %>% as.integer()) %>%
    tidyr::separate(col = id, into = c("arm", "within_arm_index"), sep = "_", remove = FALSE, convert = TRUE) %>%
    dplyr::relocate(chr, .after = arm) 
  
  tibble_wide$arm <- factor(tibble_wide$arm, levels = tibble_wide$arm %>% unique()) 
  tile_size_levels <- paste0("n_100kb_bins_", tibble_wide %>% dplyr::pull(tile_size_label) %>% sort() %>% unique() )  
  tibble_wide$tile_size <-  factor(tibble_wide$tile_size, levels = tile_size_levels)
  
  # arrange columns in specific order
  tibble_wide <- tibble_wide %>%
    dplyr::arrange(tile_size, id, within_arm_index, frag_len) %>%
    # careful
    dplyr::group_by(tile_size, frag_len) %>% 
    dplyr::mutate(id_rename = row_number() ) %>% 
    dplyr::relocate(id_rename, .after = id) %>%
    dplyr::relocate(tile_size_label, .after = tile_size) %>%
    dplyr::ungroup()
  
    return(tibble_wide)
}
