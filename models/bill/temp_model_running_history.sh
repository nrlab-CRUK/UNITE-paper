sbatch --time=1-0 --mem=64G --partition="cuda" --gres=gpu:1 -J "cnn" --wrap="singularity exec --nv --bind /scratchc/nrlab/wang04/ --bind ~ --bind /home/nrlab/wang04/ulyses/ /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python /home/nrlab/wang04/ulyses/models/bill/haichao_efficientnet_model.py --scale_pixel_to_255 --remove_problem_bin --tf_label 0.1_0.2  "

sbatch --time=1-0 --mem=64G --partition="cuda" --gres=gpu:1 -J "cnn" --wrap="singularity exec --nv --bind /scratchc/nrlab/wang04/ --bind ~ --bind /home/nrlab/wang04/ulyses/ /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python /home/nrlab/wang04/ulyses/models/bill/haichao_efficientnet_model.py  --scale_pixel_to_255 --remove_problem_bin  --tf_label 0.2_1  "



# here I changed to a higher learing rate
sbatch -J 'rate'  --output="stdout_efficientNets_tf_0_1.%j.log" --error="stdout_efficientNets_tf_0_1.%j.log"  --time=1-0 --mem=64G --partition="cuda" --gres=gpu:1 -J "cnn" --wrap="singularity exec --nv --bind /scratchc/nrlab/wang04/ --bind ~ --bind /home/nrlab/wang04/ulyses/ /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python /home/nrlab/wang04/ulyses/models/bill/haichao_efficientnet_model.py  --scale_pixel_to_255 --remove_problem_bin --model_learning_rate 0.009 --model_batch_size 60  --tf_label 0_1 --input_data_path /scratchc/nrlab/wang04/ulyses/cnn_model/cnn_strat/finalize/tp1/   "






sbatch --time=1-0 --mem=64G --partition="cuda" --gres=gpu:1 -J "cnn" --wrap="singularity exec --nv --bind /scratchc/nrlab/wang04/ --bind ~ --bind /home/nrlab/wang04/ulyses/ /scratchc/nrlab/wang04/docker/tensorflow_gpu_cuda.sif python /home/nrlab/wang04/ulyses/models/bill/haichao_efficientnet_model.py --output_path /scratchc/nrlab/wang04/ulyses_results_iteration/iteration_1_for_manuscript/cnn_model/plasma_tp1_0.1x_model_C2/output  "

