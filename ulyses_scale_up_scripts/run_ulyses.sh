#!/bin/bash
source ~/.bashrc
conda activate R4_1

# the $1 is downsampled_dir, default is pwd
downsampled_dir=$1
if [ -z "$downsampled_dir" ]; then
    downsampled_dir=$(pwd)
fi

# make a parameter called input_file_type, infer from the downsampled_dir
input_file_type="bam"
if [[ $downsampled_dir == *"finale"* ]] || [[ $downsampled_dir == *"Finale"* ]]; then
    input_file_type="rds"
fi


Rscript /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/sbatch_ulyses.R \
  --bamfile_path $downsampled_dir \
  --input_file_type ${input_file_type} \
  --replace_output FALSE \
  --block_size 2000000 \
  --slurm_mem "16G" \
  --slurm_time "0-3" \
  --slurm_job_array_task_limit 100 \
  --clear_cache FALSE

