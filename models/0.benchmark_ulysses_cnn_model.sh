# !/bin/bash
# iteration_no=8
# cnn_model="C5"  # Define your cnn_model variable here
# working_dir_base="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/"
# working_dir="${working_dir_base}cnn_model/cnn${iteration_no}/plasma_tp1_0.1x_model_${cnn_model}/"
# create the working directory if it does not exist
# mkdir -p "${working_dir}"
# cd "${working_dir}"
# mae_file="/scratchc/nrlab/wang04/build_mae/colData_batch1_2.rds_author_all_MAE.rds"
# model_epochs=40
# model_validation_split=0.30
# model_learning_rate=0.0005
# sbatch_params="--time=2-1:00:00 --mem=200G --partition=rocm --gres=gpu:2"

################################################################################
# prepare the tensor
################################################################################
#!/bin/bash
iteration_no=14
working_dir_base="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/"
mae_file="/scratchc/nrlab/wang04/build_mae/colData_batch1_2.rds_author_all_MAE.rds"
common_params="--mae_file ${mae_file} --sample_type plasma --ichorcna_tf_cutoff_lower 0 0.03 0.1 0 --ichorcna_tf_cutoff_upper 0.03 0.1 1 1 --timepoints 1"

for cnn_model in C1 C2 C3 C4 C5; do
    working_dir="${working_dir_base}cnn_model/cnn${iteration_no}/plasma_tp1_0.1x_model_${cnn_model}/"
    
    # Create the working directory if it does not exist
    mkdir -p "${working_dir}"
    cd "${working_dir}"
    
    case "$cnn_model" in
        "C1")
            layers="n_isize"
            mem="300G"
            ;;
        "C2")
            layers="n_isize n_motif_smono1_C n_motif_smono1_T"
            mem="300G"
            ;;
        "C3")
            layers="n_isize n_motif_umono1_C n_motif_umono1_T n_motif_smono1_C n_motif_smono1_T"
            mem="480G"
            ;;
        "C3_mismatch")
            layers="n_isize n_motif_umono1_C n_motif_umono1_T n_motif_smono1_C n_motif_smono1_T n_pos_nm n_neg_nm"
            mem="480G"
            authors="--authors A.S.\ et\ al A.Z.\ et\ al P.U.\ et\ al F.M.\ et\ al K.H.\ et\ al P.M.\ et\ al P.P.\ et\ al"
            ;;
        "C4")
            layers="n_isize n_motif_smono1_C n_motif_smono1_T n_motif_smono2_C n_motif_smono2_T n_motif_smono3_C n_motif_smono3_T"
            mem="480G"
            ;;
        "C5")
            layers="n_isize n_motif_umono2_C n_motif_umono2_T n_motif_umono1_C n_motif_umono1_T n_motif_smono1_C n_motif_smono1_T n_motif_smono2_C n_motif_smono2_T"
            mem="480G"
            ;;
        "C5_mismatch")
            layers="n_isize n_motif_umono3_C n_motif_umono3_T n_motif_umono2_C n_motif_umono2_T n_motif_umono1_C n_motif_umono1_T n_motif_smono1_C n_motif_smono1_T n_motif_smono2_C n_motif_smono2_T n_motif_smono3_C n_motif_smono3_T n_pos_nm n_neg_nm"
            mem="480G"
            authors="--authors A.S.\ et\ al A.Z.\ et\ al P.U.\ et\ al F.M.\ et\ al K.H.\ et\ al P.M.\ et\ al P.P.\ et\ al"
            ;;
        *)
            echo "Invalid cnn_model: $cnn_model"
            exit 1
            ;;
    esac
    
    # Submit the job using sbatch
    sbatch --time=1-0 --mem=$mem --partition="epyc" -J "${cnn_model}_tensor" \
        --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R \
        ${common_params} --outdir ${working_dir} --layers ${layers} ${authors}"
    
    # Unset authors to avoid errors in next iterations
    unset authors
done


################################################################################
# nested cv
################################################################################
paper_iteration="iteration_3_merge_below_3"
i="cnn14"
random_search_iter=10
top_layers="original"
seed=1
#which_model='resnet18'
#which_model='resnet34'
# which_model='tf_efficientnetv2_s.in1k'
# which_models=("resnet18" "resnet34" "tf_efficientnetv2_s.in1k" "plain_cnn" "vit_base")
which_models=("resnet18" "resnet34" "tf_efficientnetv2_s.in1k")
# which_models=("resnet18")
model_names=("C1" "C2" "C3" "C4" "C5")
tf_labels=("0_0.03" "0.03_0.1" "0.1_1" "0_1")
cuda_sbatch_params="--time=3-1:00:00 --mem=48G --partition=cuda --gres=gpu:1"

for which_model in "${which_models[@]}"; do
    for model_name in "${model_names[@]}"; do
        for tf_label in "${tf_labels[@]}"; do
            input_data_path="/scratchc/nrlab/wang04/ulyses_results_iteration/${paper_iteration}/cnn_model/$i/plasma_tp1_0.1x_model_${model_name}/"
            command="source ~/.bashrc; singularity exec --nv \
            --bind /scratchc/nrlab/wang04/  \
            --bind ~  \
            --bind /home/nrlab  \
            --bind /home/nrlab/wang04/ulyses/    \
            --bind /scratchc/nrlab/wang04/ulyses  \
            --bind /scratchc/nrlab/wang04/ulyses_results_iteration   \
            /scratchc/nrlab/wang04/docker/pytorch_cuda.sif python \
            /home/nrlab/wang04/ulyses/models/cnn_binary_timm_skorch_nested_cv.py \
            --tf_label ${tf_label} \
            --input_data_path ${input_data_path} \
            --remove_problem_bin \
            --model_batch_size 50  \
            --model_learning_rate 0.005  \
            --train_from_scratch \
            --model_epochs 30 \
            --hyperparameter_search_strategy random \
            --random_search_iter ${random_search_iter} \
            --inner_cv_fold 3 \
            --top_layers ${top_layers} \
            --model_name ${which_model} \
            --seed ${seed}"
            sbatch --wrap="${command}" ${cuda_sbatch_params} -J cv.${model_name}.${tf_label}.${which_model}  --output=${input_data_path}/cnn_${model_name}.${tf_label}.${which_model}_%j.out --error=${input_data_path}/cnn_${model_name}.${tf_label}.${which_model}.err
        done
    done
done

################################################################################
# train and test on independent testing data
################################################################################

paper_iteration="iteration_3_merge_below_3"
i="cnn14"
random_search_iter=10
top_layers="original"
seed=1
# which_models=("resnet18" "resnet34" "tf_efficientnetv2_s.in1k" "plain_cnn" "vit_base")
which_models=("resnet18" "resnet34" "tf_efficientnetv2_s.in1k")
# which_models=("resnet18")
model_names=("C1" "C2" "C3" "C4" "C5")
tf_labels=("0_0.03" "0.03_0.1" "0.1_1" "0_1")
cuda_sbatch_params="--time=3-1:00:00 --mem=48G --partition=cuda --gres=gpu:1"

for which_model in "${which_models[@]}"; do
    for model_name in "${model_names[@]}"; do
        for tf_label in "${tf_labels[@]}"; do
            input_data_path="/scratchc/nrlab/wang04/ulyses_results_iteration/${paper_iteration}/cnn_model/$i/plasma_tp1_0.1x_model_${model_name}/"
            command="source ~/.bashrc; singularity exec --nv \
            --bind /scratchc/nrlab/wang04/  \
            --bind ~  \
            --bind /home/nrlab  \
            --bind /home/nrlab/wang04/ulyses/    \
            --bind /scratchc/nrlab/wang04/ulyses  \
            --bind /scratchc/nrlab/wang04/ulyses_results_iteration   \
            /scratchc/nrlab/wang04/docker/pytorch_cuda.sif python \
            /home/nrlab/wang04/ulyses/models/cnn_binary_timm_skorch_train_and_test.py \
            --tf_label ${tf_label} \
            --input_data_path ${input_data_path} \
            --remove_problem_bin \
            --model_batch_size 50  \
            --model_learning_rate 0.005  \
            --train_from_scratch \
            --model_epochs 30 \
            --hyperparameter_search_strategy random \
            --random_search_iter ${random_search_iter} \
            --inner_cv_fold 3 \
            --top_layers ${top_layers} \
            --model_name ${which_model} \
            --seed ${seed}"
            sbatch --wrap="${command}" ${cuda_sbatch_params} -J test.${model_name}.${tf_label}.${which_model}  --output=${input_data_path}/cnn_${model_name}.${tf_label}.${which_model}_%j.out --error=${input_data_path}/cnn_${model_name}.${tf_label}.${which_model}.err
        done
    done
done




################################################################################
# viz the results
################################################################################

sbatch --time=00:20:00 \
    --mem=4G \
    --partition="epyc" \
    --wrap="source ~/.bashrc; conda activate R4_3; \
            Rscript /home/nrlab/wang04/ulyses/models/cnn_binary_performance_ranking_plot.R "
