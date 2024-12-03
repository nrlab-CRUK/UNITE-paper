
bam_to_image <- function(bamfile,

	# functions
	bin_functions = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/bin_functions.r",
	frag_functions = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/frag_functions.r",
	overlap_functions = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/overlap_functions.r",
	python_functions = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/functions/Python_internal_functions.py",

	# essential resource files 
	ref_csv_file = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_packaged/resources/overlap_ref.csv",
	bins_anno_tiled_file = "/mnt/scratcha/nrlab/wang04/urine/delfi/down_sample/upgrade_fragmentim/explore_QDNAseq_blacklist/bins_anno_tiled.rds",

	# output paths
	output_path_frag = dirname(bamfile),
	output_path_frag_bin_overlap = dirname(bamfile),

	# frag annoation params
	which_genome = "hg19",
	fragQC_isize_min = 0,
	fragQC_isize_max = 2000,

	# bin-frag overlap params
	isize_from = 20,
	isize_to = 500,
	overlap_summary_jobname = "overlap_summary",

	# image embedding params
	layers = c("n_isize",
            "n_motif_s1_C",
            "n_motif_s1_A",
            "n_motif_s1_G",
            "n_motif_s1_T",
            "n_neg_nm", 
            "n_pos_nm", 
            "bin_mean_GC", 
            "bin_mean_mappability")
	) {
		################################################################
		## libs and sources
		################################################################

		library(tidyverse)
		library(magrittr)
		library(Biobase)
		library(QDNAseq)
		library(Homo.sapiens)
		library(plyranges)
		library(GenomicRanges)
		library(GenomicAlignments)
		library(Gviz)
		library(matrixStats)
		library(AnnotationHub)
		library(patchwork)
		library(gganimate)
		library(rslurm)
		library(reticulate)

		source(bin_functions)
		source(frag_functions)
		source(overlap_functions)
		reticulate::source_python(python_functions)

		################################################################
		## params QC
		################################################################

		if (missing(bamfile)) {
			stop("Please indicate which bam file!") 
		}
				
		f_valid <- c(bin_functions, 
				bamfile,
				frag_functions,
				overlap_functions,
				ref_csv_file,
				bins_anno_tiled_file,
				python_functions
				)

		if(!purrr::every(f_valid, file.exists)){
			message("File(s) missing:")
			message(f_valid[-which(file.exists(f_valid))])
			message("*******************************************")
			stop("Stopped due to missing file(s).")
		}

		if (!dir.exists(output_path_frag)) {
			stop("output_path_frag: directory doesn't exist!")
		}

		if (!dir.exists(output_path_frag_bin_overlap)) {
			stop("output_path_frag_bin_overlap: directory doesn't exist!")
		}

		if (which_genome == "hg19") {
			library(BSgenome.Hsapiens.UCSC.hg19)
			genome <- BSgenome.Hsapiens.UCSC.hg19 
		} else if (which_genome == "hg38") {
			library(BSgenome.Hsapiens.UCSC.hg38)
			genome <- BSgenome.Hsapiens.UCSC.hg38 
		}

		################################################################
		## Output filenames 
		################################################################

		# frag rds file
		frag_file_basename <- paste0(basename(bamfile), ".frag.rds")
		frag_file_path <- output_path_frag
		frag_file_fullname <- file.path(frag_file_path, frag_file_basename)

		# frag_bin_overlap rds file
		overlap_file_basename <- paste0(basename(bamfile), ".frag_bin_overlap.rds")
		overlap_file_path <- output_path_frag_bin_overlap
		overlap_file_fullname <- file.path(overlap_file_path, overlap_file_basename)

		# rslurm tmp folder
		mainDir <- output_path_frag_bin_overlap 
		subDir <-  paste0(basename(bamfile), "_rslurm_tmp_dir")
		tmpDir <- file.path(mainDir, subDir) 

		# wide tibble rds file
		tibble_wide_file_basename <- paste0(basename(bamfile), ".summary_wide.rds")
		tibble_wide_file_path <- output_path_frag_bin_overlap
		tibble_wide_file_fullname <- file.path(tibble_wide_file_path, tibble_wide_file_basename)

		# long tibble rds file
		tibble_long_file_basename <- paste0(basename(bamfile), ".summary_long.rds")
		tibble_long_file_path <- output_path_frag_bin_overlap
		tibble_long_file_fullname <- file.path(tibble_long_file_path, tibble_long_file_basename)

		# long tibble csv file
		tibble_long_csv_file_basename <- paste0(basename(bamfile), ".summary_long.csv")
		tibble_long_csv_file_path <- output_path_frag_bin_overlap
		tibble_long_csv_file_fullname <- file.path(tibble_long_csv_file_path, tibble_long_csv_file_basename)

		# tenfor file 
		tensor_file_dir <- output_path_frag_bin_overlap 
		tensor_file_basename <- paste0(tibble_long_csv_file_basename, ".tensor")
		tensor_file_fullname <- file.path(tensor_file_dir, tensor_file_basename)

		if (file.exists(tensor_file_fullname)) {
			stopmsg <- paste0("tensor file exists, skip:\n", tensor_file_fullname)
			message(stopmsg)
			stop(stopmsg, call. = FALSE)
		}

		################################################################
		## Reporting parameters 
		################################################################

		message("You are calcualting these layers: \n")
		print(layers)

		message("\nFrag length ranges in the embeddings/images: \n")
		message(paste("from", isize_from, "to", isize_to, ".\n"))

		################################################################
		## function body 
		################################################################

		#### read in resource files
		ref <- readr::read_csv(ref_csv_file, show_col_types = FALSE)
		bin <- readRDS(bins_anno_tiled_file)
		n_bin <- length(bin)

		if(n_bin == 0) {

			stop("bin annotation file lenght is 0, stop running!", call. = FALSE)
		}

		

		################################################################
		## BAM to frag object
		################################################################
		if(!file.exists(frag_file_fullname)){
			frag <- bam_to_galp(bamfile = bamfile) %>%
				remove_outward_facing_readpairs() %T>%
				report_nrow(filter_step = paste0("Remove outward-facing read pairs:")) %>%
				curate_start_and_end() %T>%
				report_nrow(filter_step = paste0("Curate start and end:")) 

			seqlengths(frag) <- seqlengths(genome)[levels(seqnames(frag))]
			genome(frag) <- seqinfo(genome)@genome %>% unique()

			frag <- frag %>%
				remove_out_of_bound_reads() %T>%
				report_nrow( filter_step = paste0("Remove out-of-bound read pairs:"))

			# annotate motif

			frag <- annotate_motif(frag, genome = genome)

			# filter by frag_len, mapq, and cigar string
			frag <- frag %T>%
			report_nrow( filter_step = paste0("Before filtering:")) %>%

			plyranges::filter(frag_len >= fragQC_isize_min & frag_len <= fragQC_isize_max) %T>%
			report_nrow( filter_step = paste0("Step 1: keep isize between ", fragQC_isize_min, " and ", fragQC_isize_max, ":")) %>%

			plyranges::filter(pos_mapq >= 0 & neg_mapq >= 0) %T>%
			report_nrow( filter_step = paste0("Step 2: remove 0 map quality frags:")) %>%

			plyranges::filter(!stringr::str_detect(pos_cigar, pattern = "[ID]")) %T>%
			report_nrow( filter_step = paste0("Step 3: remove pos stand with CIGAR containing INDEL:")) %>%
			
			plyranges::filter(!stringr::str_detect(neg_cigar, pattern = "[ID]")) %T>%
			report_nrow( filter_step = paste0("Step 4: remove neg stand with CIGAR containing INDEL:"))

			# save frag rds file

			saveRDS(object = frag, file = frag_file_fullname)
			message("Saved frag file: \n")
			message(frag_file_fullname)
		} else {
			frag <- readRDS(frag_file_fullname)
			message("Skipped bin to frag step...")
		}


		################################################################
		## Overlap between bin and frag object
		################################################################
		if(!file.exists(overlap_file_fullname)){

			message("Starting bin-frag overlapping...\n")
			frag <- frag %>%
				plyranges::filter(frag_len >= isize_from & frag_len <= isize_to)

			overlap_list <- base::lapply(bin, plyranges::join_overlap_left, y = frag) %>%
			base::lapply(as_tibble)

			# save bin-frag-overlap rds file

			saveRDS(overlap_list, overlap_file_fullname)
			message("Saved bin-frag overlap file: \n")
			message(overlap_file_fullname)
		} else {
			overlap_list <- readRDS(overlap_file_fullname)

			message("Skipped frag_bin_overlap step...")
		}


		################################################################
		## Summarize the overlap into tibble 
		################################################################

		# to avoid "objects not found" error:

		bin_functions <- bin_functions
		frag_functions <- frag_functions
		overlap_functions <- overlap_functions
		python_functions <- python_functions
		which_genome <- which_genome
		fragQC_isize_max <- fragQC_isize_max
		fragQC_isize_min <- fragQC_isize_min
		isize_from <- isize_from
		isize_to <- isize_to
		layers <- layers
		output_path_frag_bin_overlap <- output_path_frag_bin_overlap
		output_path_frag <- output_path_frag
		ref_csv_file <- ref_csv_file
		bins_anno_tiled_file <- bins_anno_tiled_file
		ref <- ref
		bin <- bin
		n_bin <- n_bin

		overlap_summary_jobname <- overlap_summary_jobname


		# utilize rslurm to parallelize the summary process
		global_objects <- ls(envir = globalenv()) %>% as.vector()
		
		rslurm_tmp_dir <- file.path(tmpDir, paste0("_rslurm_", overlap_summary_jobname)) 

		if (!dir.exists(rslurm_tmp_dir)) {
			dir.create(tmpDir)
			setwd(tmpDir)
			message("Submitting overlap summary jobs via rslurm.")
			sjob <- rslurm::slurm_map(overlap_list,
					f = each_bin_size,
					output_layer = layers, 
					ref = ref,
					isize_from = isize_from,
					isize_to = isize_to,
					jobname = overlap_summary_jobname,
					nodes = 20, 
					cpus_per_node = 1,
					slurm_options = list(time = '10-0',  mem = '50G'),
					global_objects = global_objects,
					submit = TRUE)

			message("Reading rslurm output files... \n")
			sjob_result <- get_slurm_out(sjob, outtype = 'raw', wait = TRUE)
			
			
		} else {
			message(paste(rslurm_tmp_dir, "exists!"))

			# check if the rslurm files are corrupt
			rds_files <- list.files(path = rslurm_tmp_dir, pattern = "results_\\d\\.RDS")
			n_file <- length(rds_files)

			size_file <- file.size(rds_files)

			zero_in_file_size <- 0 %in% size_file 

			if(n_file == n_bin & !zero_in_file_size){
				message("Skipping summary jobs on slurm.")
				setwd(tmpDir)
				sjob <- rslurm::slurm_map(overlap_list,
						f = each_bin_size,
						output_layer = layers, 
						ref = ref,
						isize_from = isize_from,
						isize_to = isize_to,
						jobname = overlap_summary_jobname,
						nodes = 20, 
						cpus_per_node = 1,
						slurm_options = list(time = '10-0',  mem = '50G'),
						global_objects = global_objects,
						submit = FALSE)

				message("Reading rslurm output files... \n")
				sjob_result <- get_slurm_out(sjob, outtype = 'raw', wait = FALSE)

			} else {
				message("Deleting corruped rslurm tmp dir. \n")
				unlink(tmpDir, recursive=TRUE)

				# re-run the job
				dir.create(tmpDir)
				setwd(tmpDir)
				message("Submitting overlap summary jobs via rslurm.")

				sjob <- rslurm::slurm_map(overlap_list,
						f = each_bin_size,
						output_layer = layers, 
						ref = ref,
						isize_from = isize_from,
						isize_to = isize_to,
						jobname = overlap_summary_jobname,
						nodes = 20, 
						cpus_per_node = 1,
						slurm_options = list(time = '10-0',  mem = '50G'),
						global_objects = global_objects,
						submit = TRUE)

				message("Reading rslurm output files... \n")
				sjob_result <- get_slurm_out(sjob, outtype = 'raw', wait = TRUE)
			}

			
		}


		################################################################
		## Save as wide-format dataframe 
		################################################################

		message("Making wide-format dataframe...")

		# combine results into wide tibble
		tibble_wide <-  dplyr::bind_rows(sjob_result, .id = "tile_size") %>%
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

		# save summary wide tibble rds file

		message("Saving wide summary tibble: \n")
		message(tibble_wide_file_fullname)

		saveRDS(object = tibble_wide, file = tibble_wide_file_fullname )


		################################################################
		## Save as long-format dataframe 
		################################################################
		message("Making long-format dataframe...")


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

		# save summary long tibble rds file

		message("Saving long summary tibble: \n")
		message(tibble_long_file_fullname)
		saveRDS(object = tibble_long, file = tibble_long_file_fullname )

		# save as .csv file for python usage


		message("Saving long summary tibble as CSV file: \n")
		message(tibble_long_csv_file_fullname)
		readr::write_csv(x = tibble_long, file = tibble_long_csv_file_fullname )



		################################################################
		## python: convert to ndarray and save the object
		################################################################

		# functions already sourced in the beginning, ready to use
		message("Making image file via python... ")

		df_to_tensor(csv_file = tibble_long_csv_file_fullname, 
			output_dir= tensor_file_dir)
		
		################################################################
		## Cleanup works 
		################################################################
		setwd(output_path_frag_bin_overlap)

		# delete the tmp files from rslurm
		message("Deleting rslurm tmp dir: \n")
		message(tmpDir)
		unlink(tmpDir, recursive=TRUE)

		message("rslurm tmp dir deleted!")
		message("Jobs done.")

}













