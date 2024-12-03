
# train using all samples

## select all plasma samples
```bash
mkdir -p /Users/wang04/cnn_all

sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1  --outdir /Users/wang04/cnn_all "

```


## train the model

```bash

cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py /Users/wang04/cnn_all/
#cd /Users/wang04/cnn_all/

sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 ./cnn_binary_train_using_all_data.py --wd /Users/wang04/cnn_all --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

```

## train the model using VA and CS data

```bash
wd=/Users/wang04/cnn_all_VA_CS
cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${wd}


sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

```




## tensor board


```bash

ssh -L 6006:127.0.0.1:6006 wang04@clust1-sub.cri.camres.org

cd /Users/wang04/cnn_all_timepoint1

conda activate R4_2; tensorboard --logdir=logs

```


## apply the model to pp et al 

### pp et al data 

```
mkdir -p /Users/wang04/cnn_all_PP
sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1  --outdir /Users/wang04/cnn_all_PP --authors 'P.P. et al' "

```

### pp et al tp1 data
```
mkdir -p /Users/wang04/cnn_all_PP_tp1
sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1  --outdir /Users/wang04/cnn_all_PP_tp1 --authors 'P.P. et al' --timepoints 1 --sample_type plasma "

```

### fm et al tp1 data

```
wd=/Users/wang04/cnn_all_FM_tp1
mkdir -p ${wd}
sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1  --outdir ${wd} --authors 'F.M. et al' --timepoints 1 --sample_type plasma "

```


### sc et al tp1 data

```
wd=/Users/wang04/cnn_all_SC_tp1
mkdir -p ${wd}
sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1  --outdir ${wd} --authors 'S.C. et al' --timepoints 1 --sample_type plasma "

```



### load and evaluate the model using PP et al data 

```bash

cd /Users/wang04/cnn_all_PP

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 ./cnn_binary_test_using_datasets.py --wd /Users/wang04/cnn_all_PP  --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file /Users/wang04/cnn_all/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"



wd=/Users/wang04/cnn_all_PP_tp1
cd ${wd}


sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 ./cnn_binary_test_using_datasets.py --wd ${wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file /Users/wang04/cnn_all/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"


wd=/Users/wang04/cnn_all_FM_tp1
cd ${wd}


sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file /Users/wang04/cnn_all/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"





```


# train using VA, CS, test on FMtp1


```

dir=/Users/wang04/cnn_all_VA_CS
mkdir -p ${dir}
cd ${dir}
sbatch /home/nrlab/wang04/ulyses/models/select_data.sh


# train

train_wd=/Users/wang04/cnn_all_VA_CS
cd ${train_wd}
cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}



sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 40 --model_validation_split 0.40 --model_learning_rate 0.0005 --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"


#test 
test_wd=/Users/wang04/cnn_all_FM_tp1
test_wd=/Users/wang04/cnn_all_PP_tp1

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

# test using VA+CS model
sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file /Users/wang04/cnn_all_VA_CS/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"


# test using VA+CS+FM model
sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file /Users/wang04/cnn_all/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"


```
![image](https://github.com/hw538/ulyses/assets/15274940/c49d672a-bbcb-4fff-b4a7-08661120e164)



# train using VA + FM, test on SC tp1


```
# get the data 

dir=/Users/wang04/cnn_all_VA_FM
mkdir -p ${dir}
cd ${dir}
sbatch /home/nrlab/wang04/ulyses/models/select_data.sh



# train

train_wd=/Users/wang04/cnn_all_VA_FM
cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}
cd ${train_wd}

sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 40 --model_validation_split 0.40 --model_learning_rate 0.0005 --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

# test

test_wd=/Users/wang04/cnn_all_SC_tp1
cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"



```

![image](https://github.com/hw538/ulyses/assets/15274940/4cb890b1-97e6-4af1-b6ab-cdc464dd579f)

```
# test using PP_tp1

train_wd=/Users/wang04/cnn_all_VA_FM

test_wd=/Users/wang04/cnn_all_PP_tp1

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"

```


![image](https://github.com/hw538/ulyses/assets/15274940/5f2384b5-e9e8-4014-98c4-06842b62767c)




############################################################
# train only using timepoint 1 (VA, FM, SC), test on PP_tp1
############################################################

```bash

# get the data

wd=/Users/wang04/cnn_all_timepoint1
mkdir -p ${wd}
cd ${wd}

# authors_default <- c( "V.A. et al", "F.M. et al", "S.C. et al")
sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1 --timepoints 1 --sample_type plasma  --outdir ${wd} "
==============================

# train

train_wd=/Users/wang04/cnn_all_timepoint1
cd ${train_wd}
cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}



# increase epoch and learning rate
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 30 --model_validation_split 0.40 --model_learning_rate 0.001 --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

# further increase epoch
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 40 --model_validation_split 0.40 --model_learning_rate 0.001 --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"



## simpler model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 30 --model_validation_split 0.30 --model_learning_rate 0.0008 --model_type simpler --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"


:'
x_file_path:  /Users/wang04/cnn_all_timepoint1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (1058, 464, 201, 3)
y shape:  (1058,)
z shape:  (1058,)
class weights:  {0: 1.7175324675324675, 1: 0.7053333333333334}
----------------------------------------
x_train shape: (634, 464, 201, 3)
x_val shape: (424, 464, 201, 3)
y_train shape: (634,)
y_val shape: (424,)
Class balance in train set: (array([0, 1]), array([185, 449]))
Class balance in validation set: (array([0, 1]), array([123, 301]))
Epoch 1/30

- 3s 131ms/step - loss: 0.1218 - accuracy: 0.9543 - AUROC: 0.9942 - AUPRC: 0.9977 - PPV: 0.9930 - recall_sens: 0.9421 - get_f1: 0.9645 - sen_95spe: 0.9688 - sen_98spe: 0.9465 - sen_99spe: 0.9376 - val_loss: 0.2643 - val_accuracy: 0.8797 - val_AUROC: 0.9566 - val_AUPRC: 0.9840 - val_PPV: 0.9630 - val_recall_sens: 0.8638 - val_get_f1: 0.9131 - val_sen_95spe: 0.8638 - val_sen_98spe: 0.7807 - val_sen_99spe: 0.7641
Model saved to: /Users/wang04/cnn_all_timepoint1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
'


## deeper model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 35 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type deeper --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

:'

3s 130ms/step - loss: 0.1816 - accuracy: 0.9243 - AUROC: 0.9838 - AUPRC: 0.9938 - PPV: 0.9740 - recall_sens: 0.9176 - get_f1: 0.9436 - sen_95spe: 0.9131 - sen_98spe: 0.9065 - sen_99spe: 0.8753 - val_loss: 0.3095 - val_accuracy: 0.8538 - val_AUROC: 0.9575 - val_AUPRC: 0.9839 - val_PPV: 0.9799 - val_recall_sens: 0.8106 - val_get_f1: 0.8851 - val_sen_95spe: 0.8505 - val_sen_98spe: 0.7641 - val_sen_99spe: 0.7608

'

===============================
# test using PP_tp1

train_wd=/Users/wang04/cnn_all_timepoint1
test_wd=/Users/wang04/cnn_all_PP_tp1

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"

:'
------------------------------


x_file_path:  /Users/wang04/cnn_all_PP_tp1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (148, 464, 151, 3)
y shape:  (148,)
z shape:  (148,)
model_file:  /Users/wang04/cnn_all_timepoint1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
Loading model...
Evaluating...
5/5 - 10s - loss: 0.3138 - accuracy: 0.8311 - AUROC: 0.8768 - AUPRC: 0.9789 - PPV: 0.8976 - recall_sens: 0.9048 - get_f1: 0.8992 - sen_95spe: 0.7698 - sen_98spe: 0.7063 - sen_99spe: 0.7063 - 10s/epoch - 2s/step
{'test_loss': 0.31376978754997253, 'test_accuracy': 0.8310810923576355, 'test_AUROC': 0.8768036961555481, 'test_AUPRC': 0.9788652062416077, 'test_PPV': 0.8976377844810486, 'test_recall_sens': 0.9047619104385376, 'test_get_f1': 0.8991979360580444, 'test_sen_95spe': 0.7698412537574768, 'test_sen_98spe': 0.7063491940498352, 'test_sen_99spe': 0.7063491940498352}
Results saved to:  /Users/wang04/cnn_all_PP_tp1/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  0.313770
1     test_accuracy  0.831081
2        test_AUROC  0.876804
3        test_AUPRC  0.978865
4          test_PPV  0.897638
5  test_recall_sens  0.904762
6       test_get_f1  0.899198
7    test_sen_95spe  0.769841
8    test_sen_98spe  0.706349
9    test_sen_99spe  0.706349
Done!
-------------------------------



[wang04@clust1-sub-1 cnn_all_PP_tp1]$ cat slurm-29228640.out
x_file_path:  /Users/wang04/cnn_all_PP_tp1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (148, 464, 201, 3)
y shape:  (148,)
z shape:  (148,)
model_file:  /Users/wang04/cnn_all_timepoint1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
Loading model...
Evaluating...
5/5 - 10s - loss: 0.5075 - accuracy: 0.7162 - AUROC: 0.7832 - AUPRC: 0.9599 - PPV: 0.9118 - recall_sens: 0.7381 - get_f1: 0.8195 - sen_95spe: 0.5635 - sen_98spe: 0.5159 - sen_99spe: 0.5159 - 10s/epoch - 2s/step
{'test_loss': 0.5075497627258301, 'test_accuracy': 0.7162162065505981, 'test_AUROC': 0.783189058303833, 'test_AUPRC': 0.9598623514175415, 'test_PPV': 0.9117646813392639, 'test_recall_sens': 0.738095223903656, 'test_get_f1': 0.8195122480392456, 'test_sen_95spe': 0.5634920597076416, 'test_sen_98spe': 0.5158730149269104, 'test_sen_99spe': 0.5158730149269104}
Results saved to:  /Users/wang04/cnn_all_PP_tp1/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  0.507550
1     test_accuracy  0.716216
2        test_AUROC  0.783189
3        test_AUPRC  0.959862
4          test_PPV  0.911765
5  test_recall_sens  0.738095
6       test_get_f1  0.819512
7    test_sen_95spe  0.563492
8    test_sen_98spe  0.515873
9    test_sen_99spe  0.515873
Done!

'


# test using VA

train_wd=/Users/wang04/cnn_all_timepoint1
test_wd=/Users/wang04/cnn_all_VA

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"

----------------------------
x_file_path:  /Users/wang04/cnn_all_VA/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (1099, 464, 151, 3)
y shape:  (1099,)
z shape:  (1099,)
model_file:  /Users/wang04/cnn_all_timepoint1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
Loading model...
Evaluating...
35/35 - 11s - loss: 0.1634 - accuracy: 0.9409 - AUROC: 0.0000e+00 - AUPRC: 1.0000 - PPV: 1.0000 - recall_sens: 0.9409 - get_f1: 0.9697 - sen_95spe: 0.0000e+00 - sen_98spe: 0.0000e+00 - sen_99spe: 0.0000e+00 - 11s/epoch - 322ms/step
{'test_loss': 0.16343791782855988, 'test_accuracy': 0.9408553242683411, 'test_AUROC': 0.0, 'test_AUPRC': 1.0, 'test_PPV': 1.0, 'test_recall_sens': 0.9408553242683411, 'test_get_f1': 0.9696663022041321, 'test_sen_95spe': 0.0, 'test_sen_98spe': 0.0, 'test_sen_99spe': 0.0}
Results saved to:  /Users/wang04/cnn_all_VA/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  0.163438
1     test_accuracy  0.940855
2        test_AUROC  0.000000
3        test_AUPRC  1.000000
4          test_PPV  1.000000
5  test_recall_sens  0.940855
6       test_get_f1  0.969666
7    test_sen_95spe  0.000000
8    test_sen_98spe  0.000000
9    test_sen_99spe  0.000000
Done!
-----------------------------------
# test using VA tp1
train_wd=/Users/wang04/cnn_all_timepoint1
test_wd=/Users/wang04/cnn_all_VA_tp1

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"

-------------------------------------
x_file_path:  /Users/wang04/cnn_all_VA_tp1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (323, 464, 151, 3)
y shape:  (323,)
z shape:  (323,)
model_file:  /Users/wang04/cnn_all_timepoint1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
Loading model...
Evaluating...
11/11 - 10s - loss: 0.1220 - accuracy: 0.9598 - AUROC: 0.0000e+00 - AUPRC: 1.0000 - PPV: 1.0000 - recall_sens: 0.9598 - get_f1: 0.9808 - sen_95spe: 0.0000e+00 - sen_98spe: 0.0000e+00 - sen_99spe: 0.0000e+00 - 10s/epoch - 910ms/step
{'test_loss': 0.12199240177869797, 'test_accuracy': 0.9597523212432861, 'test_AUROC': 0.0, 'test_AUPRC': 0.9999999403953552, 'test_PPV': 1.0, 'test_recall_sens': 0.9597523212432861, 'test_get_f1': 0.980817437171936, 'test_sen_95spe': 0.0, 'test_sen_98spe': 0.0, 'test_sen_99spe': 0.0}
Results saved to:  /Users/wang04/cnn_all_VA_tp1/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  0.121992
1     test_accuracy  0.959752
2        test_AUROC  0.000000
3        test_AUPRC  1.000000
4          test_PPV  1.000000
5  test_recall_sens  0.959752
6       test_get_f1  0.980817
7    test_sen_95spe  0.000000
8    test_sen_98spe  0.000000
9    test_sen_99spe  0.000000
Done!

----------------------------------------

```
###########################################################
# train using FM_PP_tp1, test on SC tp1 \[*doesn't work*\]
###########################################################
```bash


# train 
train_wd=/Users/wang04/cnn_all_FM_PP_tp1
cd $train_wd

cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}



# increase epoch and learning rate

## simpler model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 35 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type simpler --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

## deeper model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 30 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type deeper --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

# it seems deeper model works better in this case
:'2s 143ms/step - loss: 0.3186 - accuracy: 0.8526 - AUROC: 0.9344 - AUPRC: 0.9587 - PPV: 0.9030 - recall_sens: 0.8359 - get_f1: 0.8653 - sen_95spe: 0.7422 - sen_98spe: 0.6680 - sen_99spe: 0.6289 - val_loss: 0.3139 - val_accuracy: 0.8707 - val_AUROC: 0.9368 - val_AUPRC: 0.9585 - val_PPV: 0.9290 - val_recall_sens: 0.8421 - val_get_f1: 0.8903 - val_sen_95spe: 0.7895 - val_sen_98spe: 0.5673 - val_sen_99spe: 0.4561
Model saved to: /Users/wang04/cnn_all_FM_SC_tp1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5'



```

###################################################################
# train only using timepoint 1 (FM, SC), test on PP_tp1 and VA_tp1
###################################################################
it seems deeper model works better in this case
trained model should be used for fpDBS test

```bash

# get the data

wd=/Users/wang04/cnn_all_FM_SC_tp1
mkdir -p ${wd}
cd ${wd}
cp /home/nrlab/wang04/ulyses/models/scale_up_binary_tensors.R  ${wd}

# authors_default <- c("F.M. et al", "S.C. et al")
sbatch --time=1-0 --mem=300G --partition="epyc" -J "tensor" --wrap="source ~/.bashrc; conda activate R4_2; Rscript scale_up_binary_tensors.R --ichorcna_tf_cutoff_lower 0  --ichorcna_tf_cutoff_upper 1 --timepoints 1 --sample_type plasma  --outdir ${wd} "
===========================


# train 
train_wd=/Users/wang04/cnn_all_FM_SC_tp1
cd $train_wd

cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}



# increase epoch and learning rate

## simpler model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 35 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type simpler --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

## deeper model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 30 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type deeper --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy"

# it seems deeper model works better in this case
:'2s 143ms/step - loss: 0.3186 - accuracy: 0.8526 - AUROC: 0.9344 - AUPRC: 0.9587 - PPV: 0.9030 - recall_sens: 0.8359 - get_f1: 0.8653 - sen_95spe: 0.7422 - sen_98spe: 0.6680 - sen_99spe: 0.6289 - val_loss: 0.3139 - val_accuracy: 0.8707 - val_AUROC: 0.9368 - val_AUPRC: 0.9585 - val_PPV: 0.9290 - val_recall_sens: 0.8421 - val_get_f1: 0.8903 - val_sen_95spe: 0.7895 - val_sen_98spe: 0.5673 - val_sen_99spe: 0.4561
Model saved to: /Users/wang04/cnn_all_FM_SC_tp1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5'




============================

# test using VC tp1
train_wd=/Users/wang04/cnn_all_FM_SC_tp1
test_wd=/Users/wang04/cnn_all_VA_tp1

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"

-----------------------------

# test using PP tp 1

train_wd=/Users/wang04/cnn_all_FM_SC_tp1
test_wd=/Users/wang04/cnn_all_PP_tp1

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5"


x_file_path:  /Users/wang04/cnn_all_PP_tp1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (148, 464, 151, 3)
y shape:  (148,)
z shape:  (148,)
model_file:  /Users/wang04/cnn_all_timepoint1//model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
Loading model...
Evaluating...
5/5 - 10s - loss: 0.2989 - accuracy: 0.8378 - AUROC: 0.8851 - AUPRC: 0.9794 - PPV: 0.8542 - recall_sens: 0.9762 - get_f1: 0.9110 - sen_95spe: 0.6349 - sen_98spe: 0.6270 - sen_99spe: 0.6270 - 10s/epoch - 2s/step
{'test_loss': 0.29887455701828003, 'test_accuracy': 0.837837815284729, 'test_AUROC': 0.8851009607315063, 'test_AUPRC': 0.9793952703475952, 'test_PPV': 0.8541666865348816, 'test_recall_sens': 0.976190447807312, 'test_get_f1': 0.9109519720077515, 'test_sen_95spe': 0.6349206566810608, 'test_sen_98spe': 0.6269841194152832, 'test_sen_99spe': 0.6269841194152832}
Results saved to:  /Users/wang04/cnn_all_PP_tp1/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  0.298875
1     test_accuracy  0.837838
2        test_AUROC  0.885101
3        test_AUPRC  0.979395
4          test_PPV  0.854167
5  test_recall_sens  0.976190
6       test_get_f1  0.910952
7    test_sen_95spe  0.634921
8    test_sen_98spe  0.626984
9    test_sen_99spe  0.626984
Done!
---

x_file_path:  /Users/wang04/cnn_all_PP_tp1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (148, 464, 201, 3)
y shape:  (148,)
z shape:  (148,)
model_file:  /Users/wang04/cnn_all_FM_SC_tp1/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.h5
Loading model...
Evaluating...
5/5 - 10s - loss: 0.2895 - accuracy: 0.8716 - AUROC: 0.8999 - AUPRC: 0.9827 - PPV: 0.9280 - recall_sens: 0.9206 - get_f1: 0.9233 - sen_95spe: 0.7619 - sen_98spe: 0.6429 - sen_99spe: 0.6429 - 10s/epoch - 2s/step
{'test_loss': 0.2895432412624359, 'test_accuracy': 0.8716216087341309, 'test_AUROC': 0.89989173412323, 'test_AUPRC': 0.9827234745025635, 'test_PPV': 0.9279999732971191, 'test_recall_sens': 0.920634925365448, 'test_get_f1': 0.9232832193374634, 'test_sen_95spe': 0.761904776096344, 'test_sen_98spe': 0.6428571343421936, 'test_sen_99spe': 0.6428571343421936}
Results saved to:  /Users/wang04/cnn_all_PP_tp1/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  0.289543
1     test_accuracy  0.871622
2        test_AUROC  0.899892
3        test_AUPRC  0.982723
4          test_PPV  0.928000
5  test_recall_sens  0.920635
6       test_get_f1  0.923283
7    test_sen_95spe  0.761905
8    test_sen_98spe  0.642857
9    test_sen_99spe  0.642857
Done!


```

## train on TF > 0.1, test on PP et al TP1 (not very good)

```bash

/Users/wang04/cnn_strat/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.npy



# train 
train_wd=/Users/wang04/cnn_strat
cd $train_wd

cp /home/nrlab/wang04/ulyses/models/cnn_binary_train_using_all_data.py ${train_wd}



# increase epoch and learning rate

## simpler model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 30 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type simpler --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.npy"

## deeper model
sbatch --time=0-1:00:00 --mem=200G --partition="rocm" --gres=gpu:2 --wrap="spack load singularity@3.8.5;singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3  cnn_binary_train_using_all_data.py --wd ${train_wd} --model_epochs 30 --model_validation_split 0.40 --model_learning_rate 0.001 --model_type deeper --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.npy"


# test using PP tp 1

train_wd=/Users/wang04/cnn_strat
test_wd=/Users/wang04/cnn_all_PP_tp1

cp /home/nrlab/wang04/ulyses/models/cnn_binary_test_using_datasets.py ${test_wd}
cd ${test_wd}

sbatch --time=0-1:00:00 --mem=128G --partition="rocm" --gres=gpu:1 --wrap="spack load singularity@3.8.5; singularity exec /home/software/images/rocm/tensorflow-autobuilds_latest.sif python3 cnn_binary_test_using_datasets.py --wd ${test_wd} --x_file x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy --model_file ${train_wd}/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.binary_model.h5"



[wang04@clust1-sub-1 cnn_all_PP_tp1]$ cat slurm-29226773.out
x_file_path:  /Users/wang04/cnn_all_PP_tp1/x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.npy
x shape:  (148, 464, 201, 3)
y shape:  (148,)
z shape:  (148,)
model_file:  /Users/wang04/cnn_strat/model_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.binary_model.h5
Loading model...
Evaluating...
5/5 - 10s - loss: 2.7622 - accuracy: 0.4595 - AUROC: 0.8252 - AUPRC: 0.9684 - PPV: 1.0000 - recall_sens: 0.3651 - get_f1: 0.5274 - sen_95spe: 0.6111 - sen_98spe: 0.5635 - sen_99spe: 0.5635 - 10s/epoch - 2s/step
{'test_loss': 2.7621841430664062, 'test_accuracy': 0.45945945382118225, 'test_AUROC': 0.8252164721488953, 'test_AUPRC': 0.9683666825294495, 'test_PPV': 1.0, 'test_recall_sens': 0.3650793731212616, 'test_get_f1': 0.5273876190185547, 'test_sen_95spe': 0.6111111044883728, 'test_sen_98spe': 0.5634920597076416, 'test_sen_99spe': 0.5634920597076416}
Results saved to:  /Users/wang04/cnn_all_PP_tp1/performance_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_1_.binary_model.csv
             metric    result
0         test_loss  2.762184
1     test_accuracy  0.459459
2        test_AUROC  0.825216
3        test_AUPRC  0.968367
4          test_PPV  1.000000
5  test_recall_sens  0.365079
6       test_get_f1  0.527388
7    test_sen_95spe  0.611111
8    test_sen_98spe  0.563492
9    test_sen_99spe  0.563492
Done!

```



## tensor board

```bash

ssh -L 6006:127.0.0.1:6006 wang04@clust1-sub.cri.camres.org

cd /Users/wang04/cnn_all_timepoint1

conda activate R4_2; tensorboard --logdir=logs

```




# train using all TF > 5%


