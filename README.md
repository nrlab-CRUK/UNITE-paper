The data of the UNITE paper is deposited in the `raw_data` folder in this repository.
The meta data of each sample is in the `Table S12` in the manuscript. (Aslso see the meta_for_samples_to_upload.xlsx file in the repository)

The snippets/steps of running the UNITE pipeline are as follows:
# Prepare model input features/arrays

## Step 1: index bam files (when needed)
```bash
for i in *.bam; do sbatch --time=0-1:00:00 --mem=4G --partition=epyc -J index --wrap="source ~/.bashrc; conda activate projectx; samtools index ${i}"; done
```

## Step 2: calculate bam stats (e.g., number of reads, and coverage.)
```bash
for i in *.bam
do
sbatch --time=0-10 --mem=128G  -J "count" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_2; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/summariseBam.R  ${i}  .  FALSE"
```


## Step 3: downsample the bam file to target_depth (modify the variables below when needed) 
```bash
downsampled_dir="../downsampled_bam"
#the unit of target_depth is "x"
target_depth=0.1
mkdir -p ${downsampled_dir}
for i in *.bam.summary.csv
do
sbatch --time=0-20 --mem=32G  -J "downsamp" -p "epyc" --wrap="source ~/.bashrc; conda activate R4_1; Rscript --vanilla /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/downsample_gatk.R ${i}  ${target_depth}  ${downsampled_dir}  /home/nrlab/tools/anaconda3/envs/projectx/bin/gatk"
done
```

## Step 4: run UNITE
```bash
working_dir=""
/home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/run_ulyses.sh ${working_dir}
```

## Step 5: UNITE plots

In last step, for each bam file, there will be "*ulyses_image.rds" files generated.
The code below will generate corresponding plots for each rds file generated.
```bash
for i in *.ulyses_image.rds
do 
sbatch --time=0-00:15:00 \
  --mem=1G --partition=epyc \
  -J "plot" \
  --wrap="source ~/.bashrc; \
    conda activate R4_2; \
    Rscript /home/nrlab/wang04/ulyses/explore_packaged/functions/ulyses_image_to_features.R \
      --image_file $i \
      --specimen plasma \
      --isize_filter_min 100 \
      --isize_filter_max 220 "
done
```

## Step 6: Quality Control the UNITE results
This step checks if the result files for each samples are successfully generated.
Enter the directory containing the UNITE result files, then run the code below:
```bash
sbatch --time=1-0 --mem=8G -J QC --partition=epyc --wrap="bash /home/nrlab/wang04/ulyses/ulyses_scale_up_scripts/QC_ulyses_results.sh"
```

# using the best model to examine the test data

```bash
# best CNN model trained using 0-3% TF
# model was saved using this code: https://github.com/hw538/ulyses/blob/7117a73b9c5490f40c69a2550040e75cda5c8c33/models/cnn_binary_timm_skorch_train_and_test.py#L620
/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn14/plasma_tp1_0.1x_model_C4/roc_curve_resnet18_tf_0_0.03/ind_test_random/resnet18/repeat_0_fold_0_best_model_from_CV.pkl
# best CNN model trained using all TF
/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn14/plasma_tp1_0.1x_model_C4/roc_curve_resnet18_tf_0_1/ind_test_random/resnet18/repeat_0_fold_0_best_model_from_CV.pkl
```
