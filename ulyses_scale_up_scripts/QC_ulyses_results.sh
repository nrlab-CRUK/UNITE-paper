#!/bin/bash

current_dir=$(pwd)

bam_files=($current_dir/*.bam)
ulyses_image_files=($current_dir/*.ulyses_image.rds)
ulyses_obj_files=($current_dir/*.ulyses_obj.rds)
ulyses_grid_files=($current_dir/*.ulyses_grid.pdf)
ulyses_bind_files=($current_dir/*.ulyses_bind_olap_chunks.rds)
ulyses_summary_files=($current_dir/*.ulyses_layer_summary.rds)
ulyses_plot_files=($current_dir/*.ulyses_plot.pdf)
ulyses_normalization_plot_files=($current_dir/*.ulyses_plot.pdf_CNV_normalization_plot.pdf)

num_bam_files=${#bam_files[@]}
num_ulyses_image_files=${#ulyses_image_files[@]}
num_ulyses_obj_files=${#ulyses_obj_files[@]}
num_ulyses_grid_files=${#ulyses_grid_files[@]}
num_ulyses_bind_files=${#ulyses_bind_files[@]}
num_ulyses_summary_files=${#ulyses_summary_files[@]}
num_ulyses_plot_files=${#ulyses_plot_files[@]}
num_normalization_plot_files=${#ulyses_normalization_plot_files[@]}

# Create a file to store the names of failed files or clear it if it already exists
failed_file="failed_files.txt"
if [ -e "$failed_file" ]; then
    > "$failed_file" || { echo "Error: Unable to clear $failed_file." >&2; exit 1; }
else
    touch "$failed_file" || { echo "Error: Unable to create $failed_file." >&2; exit 1; }
fi

# test if the number of files for each type is the same
if [ $num_bam_files -eq $num_ulyses_image_files ] && [ $num_ulyses_image_files -eq $num_ulyses_obj_files ] && [ $num_ulyses_obj_files -eq $num_ulyses_grid_files ] && [ $num_ulyses_grid_files -eq $num_ulyses_bind_files ] && [ $num_ulyses_bind_files -eq $num_ulyses_summary_files ] && [ $num_ulyses_summary_files -eq $num_ulyses_plot_files ] && [ $num_ulyses_plot_files -eq $num_normalization_plot_files ]; then
    echo "The number of files for each type is the same."  >> failed_files.txt
    echo "Number of bam files: $num_bam_files" >> failed_files.txt
else
    echo "The number of files for each type is NOT the same.">> failed_files.txt
    echo "Number of bam files: $num_bam_files">> failed_files.txt
    echo "Number of ulyses_image files: $num_ulyses_image_files">> failed_files.txt
    echo "Number of ulyses_obj files: $num_ulyses_obj_files">> failed_files.txt
    echo "Number of ulyses_grid files: $num_ulyses_grid_files">> failed_files.txt
    echo "Number of ulyses_bind files: $num_ulyses_bind_files">> failed_files.txt
    echo "Number of ulyses_summary files: $num_ulyses_summary_files">> failed_files.txt
    echo "Number of ulyses_plot files: $num_ulyses_plot_files">> failed_files.txt
    echo "Number of normalization_plot files: $num_normalization_plot_files">> failed_files.txt
fi




# Set the path to your R4_1 environment
conda_env="R4_1"

# Source the .bashrc only if it exists
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
fi

# Activate the conda environment
if command -v conda &> /dev/null; then
    conda activate "$conda_env" || { echo "Error: Unable to activate conda environment." >&2; exit 1; }
else
    echo "Error: conda command not found. Please make sure conda is installed and in your PATH." >&2
    exit 1
fi

# Get a list of *.ulyses_image.rds files
ulyses_image_files=(*.ulyses_image.rds)

# Function to check if an RDS file is an R list with 10 elements
check_rds_file() {
    local file="$1"
    local result=$(Rscript -e "data <- readRDS('$file'); cat(length(data) == 10)")
    if [ "$result" == "TRUE" ]; then
    # do nothing, no-op
        :
    else
        echo "$file does not meet the criteria.--------------------------------"
        echo "$file" >> failed_files.txt
    fi
}



# Iterate through each *.ulyses_image.rds file and check
for file in "${ulyses_image_files[@]}"; do
    check_rds_file "$file"
done

# Deactivate the conda environment
conda deactivate || { echo "Error: Unable to deactivate conda environment." >&2; exit 1; }
