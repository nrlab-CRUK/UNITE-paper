

# Define the parameters
mae_file="/scratchc/nrlab/wang04/ulyses/data/MAE3.rds"
model_epochs=40
model_validation_split=0.30
model_learning_rate=0.0005
sbatch_params="--time=2-1:00:00 --mem=200G --partition=rocm --gres=gpu:2"

# Define the x_files
x_files=("x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_0.01_.npy"
         "x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.01_0.03_.npy"
         "x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.03_0.1_.npy"
         "x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.1_0.2_.npy"
         "x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.npy")



############################################################################################################################################
# tp1 data
############################################################################################################################################

# prepare data
# remember to change layers and timepoints

cd /scratchc/nrlab/wang04/ulyses/cnn_model/cnn_strat/finalize/tp1
cp /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R ./
cp /home/nrlab/wang04/ulyses/models/cnn_binary_simpler.py ./

sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript scale_up_binary_tensors.R   --outdir $(pwd)  --mae_file ${mae_file} --timepoints 1"

#train and report metrics

for x_file in "${x_files[@]}"; do
    command="singularity exec --bind \$(pwd) /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_simpler.py --model_epochs ${model_epochs} --model_validation_split ${model_validation_split} --model_learning_rate ${model_learning_rate} --x_file ${x_file}"
    sbatch --wrap="${command}" ${sbatch_params}
done

# viz the results

sbatch --time=00:20:00 --mem=4G --partition="epyc" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/cnn_binary_performance_viz.R "

#-------------------------------------------------------------
 # all tp1 data [0, 1]
#-------------------------------------------------------------

cd /scratchc/nrlab/wang04/ulyses/cnn_model/cnn_strat/finalize/tp1
cp /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors_0-1.R ./
cp /home/nrlab/wang04/ulyses/models/cnn_binary_simpler.py ./

# get the data
sbatch --time=1-0 --mem=300G --partition="epyc" -J "0-1 tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript scale_up_binary_tensors_0-1.R   --outdir $(pwd)  --mae_file ${mae_file} --timepoints 1"

# run the model
# don't forget to change the the timepoint to 1
all_file="x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"
command="singularity exec --bind \$(pwd) /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_simpler.py --model_epochs ${model_epochs} --model_validation_split ${model_validation_split} --model_learning_rate ${model_learning_rate} --x_file ${all_file}"
sbatch --wrap="${command}" ${sbatch_params}

# viz the results

sbatch --time=00:20:00 --mem=4G --partition="epyc" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/cnn_binary_performance_viz.R "

############################################################################################################################################
# all tp data
############################################################################################################################################

# prepare data

cd /scratchc/nrlab/wang04/ulyses/cnn_model/cnn_strat/finalize/tp_all
cp /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R ./
cp /home/nrlab/wang04/ulyses/models/cnn_binary_simpler.py ./

sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript scale_up_binary_tensors.R   --outdir $(pwd)  --mae_file ${mae_file} "

#train and report metrics

for x_file in "${x_files[@]}"; do
    command="singularity exec --bind \$(pwd) /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_simpler.py --model_epochs ${model_epochs} --model_validation_split ${model_validation_split} --model_learning_rate ${model_learning_rate} --x_file ${x_file}"
    sbatch --wrap="${command}" ${sbatch_params}
done
# viz the results

sbatch --time=00:20:00 --mem=4G --partition="epyc" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/cnn_binary_performance_viz.R "
