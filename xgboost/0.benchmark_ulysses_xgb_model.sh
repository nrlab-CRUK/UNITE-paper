iteration_no=11
# Define the working directory
wd_base="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model"
wd="${wd_base}/xgb${iteration_no}" 
# create the working directory if it does not exist
if [ ! -d $wd ]; then
  mkdir -p $wd
fi

cd $wd

############################################################################################################################################
#  make the input data
############################################################################################################################################

tf_labels=("all"
  "0-0.03_pairwise"
  "0.03-0.1_pairwise"
  "0.1-1_pairwise")

for tf_label in "${tf_labels[@]}"; do
  sbatch --time=0-10 --mem=64G --partition="epyc" -J ${tf_label} --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/xgboost/xgboost_prepare_data.R --tf_label ${tf_label} --outdir ${wd}"
done



########################################################################################################################
# train the model
########################################################################################################################
# Use find to get a list of files and store them in an array
# cuda_sbatch_params="--time=1-1:00:00 --mem=48G --partition=cuda --gres=gpu:1 --output=${wd}/xgboost_train_%j.out --error=${wd}/xgboost_train_%j.err"

# IFS=$'\n' read -r -d '' -a files < <(find "${wd}" -type f -name 'xgboost_input.feat*.csv' && printf '\0')

# # Ensure the files array is not empty
# if [ ${#files[@]} -eq 0 ]; then
#     echo "Error: No files found."
#     exit 1
# fi
# formal run, buffered output

# for file in "${files[@]}"; do
#     bn=$(basename ${file})  # Get the basename
#     filename_without_extension=${bn%.*} 
#     command="source ~/.bashrc;singularity exec --nv \
#         --bind /scratchc/nrlab/wang04/ \
#         --bind ~ \
#         --bind /home/nrlab \
#         --bind /home/nrlab/wang04/ulyses/ \
#         --bind /scratchc/nrlab/wang04/ulyses \
#         --bind /scratchc/nrlab/wang04/ulyses_results_iteration \
#         /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python  \
#         /home/nrlab/wang04/ulyses/xgboost/xgboost_binary_function_nested_cv.py \
#         --xgb_training_device cuda \
#         --input_csv ${file}"
#     sbatch --wrap="${command}" ${cuda_sbatch_params} -J ${filename_without_extension}
# done



# run nested cv on cpu with random search for hyperparameter tuning

IFS=$'\n' read -r -d '' -a files < <(find "${wd}" -type f -name 'xgboost_input.feat*.csv' && printf '\0')

# Ensure the files array is not empty
if [ ${#files[@]} -eq 0 ]; then
    echo "Error: No files found."
    exit 1
fi

########################################################################################################################
# nested cv
########################################################################################################################

cuda_sbatch_params="--time=1-1:00:00 --mem=48G --partition=cuda --gres=gpu:1 --output=${wd}/xgboost_nested_cv_%j.out --error=${wd}/xgboost_nested_cv_%j.err"
epyc_sbatch_params="--time=3-10:00:00 --mem=48G --partition=epyc --output=${wd}/xgboost_nested_cv_%j.out --error=${wd}/xgboost_nested_cv_%j.err"
for file in "${files[@]}"; do
    bn=$(basename ${file})  # Get the basename
    filename_without_extension=${bn%.*} 
    command="source ~/.bashrc;singularity exec \
        --bind /scratchc/nrlab/wang04/ \
        --bind ~ \
        --bind /home/nrlab \
        --bind /home/nrlab/wang04/ulyses/ \
        --bind /scratchc/nrlab/wang04/ulyses \
        --bind /scratchc/nrlab/wang04/ulyses_results_iteration \
        /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python  \
        /home/nrlab/wang04/ulyses/xgboost/xgboost_binary_function_nested_cv.py \
        --xgb_training_device cpu \
        --hyperparameter_search_strategy random \
        --random_search_iter 10 \
        --input_csv ${file}"
    sbatch --wrap="${command}" ${epyc_sbatch_params} -J ${filename_without_extension}
done


########################################################################################################################
# run the testing on independent data on cpu 
########################################################################################################################

cuda_sbatch_params="--time=1-1:00:00 --mem=48G --partition=cuda --gres=gpu:1 --output=${wd}/xgboost_indtest_%j.out --error=${wd}/xgboost_indtest_%j.err"
epyc_sbatch_params="--time=1-1:00:00 --mem=48G --partition=epyc --output=${wd}/xgboost_indtest_%j.out --error=${wd}/xgboost_indtest_%j.err"
for file in "${files[@]}"; do
    bn=$(basename ${file})  # Get the basename
    filename_without_extension=${bn%.*} 
    command="source ~/.bashrc;singularity exec \
        --bind /scratchc/nrlab/wang04/ \
        --bind ~ \
        --bind /home/nrlab \
        --bind /home/nrlab/wang04/ulyses/ \
        --bind /scratchc/nrlab/wang04/ulyses \
        --bind /scratchc/nrlab/wang04/ulyses_results_iteration \
        /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python \
         /home/nrlab/wang04/ulyses/xgboost/xgboost_train_and_test_on_exam_data.py \
        --xgb_training_device cpu \
        --input_csv ${file}"
    sbatch --wrap="${command}" ${epyc_sbatch_params} -J ${filename_without_extension}
done




