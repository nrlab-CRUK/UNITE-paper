
plot_image <- function(
	# IO params 
	file,
	output_path = dirname(file),
	# filtering params
	tile_size_filter_label = c(50),
	frag_len_min = 0,
	frag_len_max = 500,
	# ggplot2 params
	plot_format = "pdf",
	x_axis_breaks = c(80, 100, 146, 166, 200, 300, 400, 500),
	x_axis_labels = c("80", "100", "146", "166", "200", "300", "400", "500"),
	plot_width = 20,
	plot_height = 20,
	plot_unit = "in", ...){ 

		################################################################
		## libs and sources
		################################################################

		library(tidyverse)
		library(magrittr)
		library(patchwork)
		library(gganimate)

		################################################################
		## params QC
		################################################################
		# check if file exist 
		if (!file.exists(file)) {
			stopmsg <- paste0("Input file doesn't exist: \n", file)
			stop(stopmsg, call. = FALSE)
		} else {
			message("Reading input file...")
			tb <- readr::read_csv(file, show_col_types = FALSE)
		}

		# check if ggplot2 params legit

		if(length(x_axis_breaks) != length(x_axis_labels)){

			stopmsg <- paste0("'x_axis_breaks'", " and ", "'x_axis_labels'", " have different number of elements!")
			stop(stopmsg, call. = FALSE)
		}

		################################################################
		## Generating plots  
		################################################################

		message("Filtering dataframe...")

		tb_filtered <- tb %>%
			dplyr::filter(tile_size_label %in% all_of(!!tile_size_filter_label)) %>%
			dplyr::filter(frag_len >= !!frag_len_min) %>%
			dplyr::filter(frag_len <= !!frag_len_max)

		message("Plotting...")

		p1 <- ggplot(tb_filtered, aes(frag_len, id_rename)) + 
			geom_raster(aes(fill = pixel)) +
			theme_classic() +
			facet_grid(tile_size_label ~ channel, scales = "free", space = "free_y") +
			scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels, expand = c(0, 0)) +
			scale_y_discrete(expand = c(0, 0)) +
			theme(axis.text.y=element_blank(),
				axis.ticks.y=element_blank(),
				axis.text.x = element_text(angle = 45, hjust=1, size = 6))

		################################################################
		## Saving plot files  
		################################################################
		
		# long tibble csv file
		tile_channel_grid_plot <- paste0(basename(file), ".grid_plot.", plot_format)
		tile_channel_grid_plot_fullname <- file.path(output_path, tile_channel_grid_plot)

		ggsave(plot = p1, 
			filename = tile_channel_grid_plot_fullname, 
			width = plot_width, 
			height = plot_height, 
			unit = plot_unit)

		message(paste("Saved plot file: \n", tile_channel_grid_plot_fullname))

}
