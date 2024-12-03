

if [ "$cnn_model" = "C1" ]; then
    working_dir="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_1_for_manuscript/cnn_model/plasma_tp1_0.1x_model_C1"
    mae_file="/scratchc/nrlab/wang04/build_mae/archive/colData_batch1_2.rds_author_all_MAE.rds"
    model_epochs=40
    model_validation_split=0.30
    model_learning_rate=0.0005
    sbatch_params="--time=2-1:00:00 --mem=200G --partition=rocm --gres=gpu:2"

    # Change working directory
    cd "${working_dir}"

    # Submit the job using sbatch
    sbatch --time=1-0 --mem=300G --partition="epyc" -J "C1_tensor" \
    --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R   \
        --outdir ${working_dir} \
        --mae_file ${mae_file} \
        --sample_type plasma \
        --layers n_isize \
        --ichorcna_tf_cutoff_lower 0  \
        --ichorcna_tf_cutoff_upper 0.01 \
        --timepoints 1"
fi


if [ "$cnn_model" = "C2" ]; then
    working_dir="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_1_for_manuscript/cnn_model/plasma_tp1_0.1x_model_C2"
    mae_file="/scratchc/nrlab/wang04/build_mae/archive/colData_batch1_2.rds_author_all_MAE.rds"
    model_epochs=40
    model_validation_split=0.30
    model_learning_rate=0.0005
    sbatch_params="--time=2-1:00:00 --mem=200G --partition=rocm --gres=gpu:2"

    # Change working directory
    cd "${working_dir}"

    # Submit the job using sbatch
    sbatch --time=1-0 --mem=300G --partition="epyc" -J "C2_tensor" \
    --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R   \
        --outdir ${working_dir} \
        --mae_file ${mae_file} \
        --sample_type plasma \
        --layers n_isize n_motif_smono1_C n_motif_smono1_T \
        --ichorcna_tf_cutoff_lower 0  \
        --ichorcna_tf_cutoff_upper 0.01 \
        --timepoints 1"
fi



if [ "$cnn_model" = "C3" ]; then
    working_dir="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_1_for_manuscript/cnn_model/plasma_tp1_0.1x_model_C3"
    mae_file="/scratchc/nrlab/wang04/build_mae/archive/colData_batch1_2.rds_author_all_MAE.rds"

    # Change working directory
    cd "${working_dir}"

    # Submit the job using sbatch
    sbatch --time=2-0 --mem=320G --partition="epyc" -J "C3_tensor" \
    --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R   \
        --outdir ${working_dir} \
        --mae_file ${mae_file} \
        --sample_type plasma \
        --layers n_isize n_motif_smono1_C n_motif_smono1_T n_motif_smono2_C n_motif_smono2_T n_motif_smono3_C n_motif_smono3_T \
        --ichorcna_tf_cutoff_lower 0 \
        --ichorcna_tf_cutoff_upper 0.01 \
        --timepoints 1"
fi

if [ "$cnn_model" = "C4" ]; then
    working_dir="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_1_for_manuscript/cnn_model/plasma_tp1_0.1x_model_C4"
    mae_file="/scratchc/nrlab/wang04/build_mae/archive/colData_batch1_2.rds_author_all_MAE.rds"
    # Change working directory
    cd "${working_dir}"

    # Submit the job using sbatch
    sbatch --time=2-0 --mem=400G --partition="epyc" -J "C4_tensor" \
    --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R   \
        --outdir ${working_dir} \
        --mae_file ${mae_file} \
        --sample_type plasma \
        --layers n_isize n_motif_umono3_C n_motif_umono3_T n_motif_umono2_C n_motif_umono2_T n_motif_umono1_C n_motif_umono1_T n_motif_smono1_C n_motif_smono1_T n_motif_smono2_C n_motif_smono2_T n_motif_smono3_C n_motif_smono3_T \
        --ichorcna_tf_cutoff_lower 0 0  \
        --ichorcna_tf_cutoff_upper 1 0.01 \
        --timepoints 1"
fi



if [ "$cnn_model" = "C5" ]; then
    working_dir="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_1_for_manuscript/cnn_model/plasma_tp1_0.1x_model_C5"
    # create the working directory if it does not exist
    mkdir -p "${working_dir}"
    mae_file="/scratchc/nrlab/wang04/build_mae/archive/colData_batch1_2.rds_author_all_MAE.rds"

    # Change working directory
    cd "${working_dir}"

    # Submit the job using sbatch
    sbatch --time=2-0 --mem=400G --partition="epyc" -J "C5_tensor" \
    --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R   \
        --outdir ${working_dir} \
        --mae_file ${mae_file} \
        --sample_type plasma \
        --layers n_isize n_motif_umono3_C n_motif_umono3_T n_motif_umono2_C n_motif_umono2_T n_motif_umono1_C n_motif_umono1_T n_motif_smono1_C n_motif_smono1_T n_motif_smono2_C n_motif_smono2_T n_motif_smono3_C n_motif_smono3_T n_pos_nm n_neg_nm \
        --ichorcna_tf_cutoff_lower  0  \
        --ichorcna_tf_cutoff_upper  0.01 \
        --authors P.P.\ et\ al F.M.\ et\ al A.Z.\ et\ al K.H.\ et\ al A.S.\ et\ al P.M.\ et\ al E.H.\ et\ al \
        --timepoints 1"
fi

