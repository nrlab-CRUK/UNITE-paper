
from sklearn.metrics import f1_score
import keras
from keras import layers
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.utils import class_weight
import sklearn.model_selection
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_curve, precision_score
from sklearn.metrics import confusion_matrix
from collections import defaultdict
import os
import argparse
import random
from sklearn.datasets import make_classification
from collections import defaultdict
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import classification_report

parser = argparse.ArgumentParser(description='Script parameters')
parser.add_argument('--input_data_path', type=str,
                    default='/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn7/plasma_tp1_0.1x_model_C2/', help='Path to input data')
parser.add_argument('--tf_label', type=str,
                    default='0_0.03', help='TF label')
parser.add_argument('--model_name', type=str,
                    default='efficient_net_s', help='Model name')
parser.add_argument('--remove_problem_bin',
                    action='store_true', help='Remove problem bin')
# parser.add_argument('--output_path', type=str,
                    # default='/home/nrlab/wang04/ulyses/models/bill/', help='Path to output files')

# parser.add_argument('--scale_pixel_to_255',
                    # action='store_true', help='Scale pixel values to 255')
parser.add_argument('--num_classes', type=int,
                    default=1, help='Number of classes')
parser.add_argument('--model_learning_rate', type=float,
                    default=0.004, help='Learning rate for the model')
parser.add_argument('--model_finetune_learning_rate', type=float,
                    default=1e-5, help='Learning rate for fine-tuning the model')
parser.add_argument('--num_repeats', type=int,
                    default=10, help='Number of repeats')
parser.add_argument('--num_folds', type=int, default=5,
                    help='Number of folds for cross-validation')
parser.add_argument('--train_val_split_fold', type=int,
                    default=3, help='Fold number for train-validation split')
parser.add_argument('--model_epochs', type=int, default=40,
                    help='Number of epochs for training the model')
parser.add_argument('--model_batch_size', type=int, default=30,
                    help='Batch size for training the model')
parser.add_argument('--model_fine_tune_epochs', type=int,
                    default=4, help='Number of epochs for fine-tuning the model')
parser.add_argument('--set_sample_weight_for_testing',
                    action='store_true', help='Set sample weight for testing')

# add a parameter called "seed" to the parser
parser.add_argument('--seed', type=int, default=0, help='Random seed')

args = parser.parse_args()

input_data_path = args.input_data_path
tf_label = args.tf_label
model_name = args.model_name
remove_problem_bin = args.remove_problem_bin
# output_path = args.output_path
# scale_pixel_to_255 = args.scale_pixel_to_255
num_classes = args.num_classes
model_learning_rate = args.model_learning_rate
model_finetune_learning_rate = args.model_finetune_learning_rate
num_repeats = args.num_repeats
num_folds = args.num_folds
train_val_split_fold = args.train_val_split_fold
model_epochs = args.model_epochs
model_batch_size = args.model_batch_size
model_fine_tune_epochs = args.model_fine_tune_epochs
if_set_sample_weight_for_testing = args.set_sample_weight_for_testing
seed = args.seed

output_path = os.path.dirname(
    input_data_path) + '/'

results = {}
results_container = defaultdict(list)

x_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'x_clean_n') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
y_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'y_clean_n') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
z_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'z_clean_n') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
bam_id_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'bam_id_clean_n') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]

x_file_path = os.path.join(input_data_path, x_file_path_file)
y_file_path = os.path.join(input_data_path, y_file_path_file)
z_file_path = os.path.join(input_data_path, z_file_path_file)
bam_id_file_path = os.path.join(input_data_path, bam_id_file_path_file)


# stop the program if x_file_path, y_file_path, z_file_path are not found
if not x_file_path or not y_file_path or not z_file_path:
    raise ValueError(
        'x_file_path, y_file_path, z_file_path are not found in the input_data_path')

# output
performance_csv_path = output_path + \
    'performance_' + model_name + '_tf_' + tf_label + '.csv'

parameters_csv_path = output_path + \
    'parameters_' + model_name + '_tf_' + tf_label + '.csv'
# roc_curve_csv_path as the same as parent dir of performance_csv_path but with different subfolder name called roc_curve
roc_curve_csv_path = output_path + 'roc_curve_' + \
    model_name + '_tf_' + tf_label + '/'

# create folder under roc_curve_csv_path called ind_test
roc_curve_csv_path_ind_test = roc_curve_csv_path + 'ind_test/'

# create the roc_curve_csv_path if not exist
os.makedirs(roc_curve_csv_path, exist_ok=True)
os.makedirs(roc_curve_csv_path_ind_test, exist_ok=True)


# read in indtest data

# indtest files
x_indtest_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'x_clean_independent_test') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
y_indtest_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'y_clean_independent_test') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
z_indtest_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'z_clean_independent_test') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
bam_id_indtest_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'bam_id_clean_independent_test') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]

x_indtest_file_path = os.path.join(input_data_path, x_indtest_file_path_file)
y_indtest_file_path = os.path.join(input_data_path, y_indtest_file_path_file)
z_indtest_file_path = os.path.join(input_data_path, z_indtest_file_path_file)
bam_id_indtest_file_path = os.path.join(
    input_data_path, bam_id_indtest_file_path_file)


# load the data
x = np.load(x_file_path, allow_pickle=True)
y = np.load(y_file_path, allow_pickle=True)
z = np.load(z_file_path, allow_pickle=True)
bam_id = np.load(bam_id_file_path, allow_pickle=True)

# load independent test data
x_indtest = np.load(x_indtest_file_path, allow_pickle=True)
y_indtest = np.load(y_indtest_file_path, allow_pickle=True)
z_indtest = np.load(z_indtest_file_path, allow_pickle=True)
bam_id_indtest = np.load(bam_id_indtest_file_path, allow_pickle=True)


# compute class weights
model_class_weight = dict(
    zip(np.unique(y),
        class_weight.compute_class_weight(class_weight='balanced',
                                          classes=np.unique(np.array(y)), y=y)))


def sensitivity_at_specificity(y_true, y_pred, specificity, sample_weight=None):
    fpr, tpr, thresholds = roc_curve(
        y_true, y_pred, pos_label=1, sample_weight=sample_weight)
    idx = np.argmin(np.abs(fpr - (1 - specificity)))
    return tpr[idx], thresholds[idx]


def ppv_at_specificity(y_true, y_pred, threshold, sample_weight=None):
    y_pred_above_threshold = y_pred > threshold
    y_pred_above_threshold = y_pred_above_threshold.astype(int)
    return precision_score(y_true, y_pred_above_threshold, sample_weight=sample_weight)

# write a function for calculating the pred accuracy at a given specificity


def accuracy_at_specificity(y_true, y_pred, specificity, sample_weight=None):
    # Calculate the threshold corresponding to the specified specificity
    fpr, _, thresholds = roc_curve(
        y_true, y_pred, pos_label=1, sample_weight=sample_weight)
    idx = np.argmin(np.abs(fpr - (1 - specificity)))
    threshold = thresholds[idx]
    # Calculate confusion matrix
    tn, fp, fn, tp = confusion_matrix(
        y_true, y_pred > threshold, sample_weight=sample_weight).ravel()
    # Calculate specificity and sensitivity
    tnr = tn / (tn + fp)
    sensitivity = tp / (tp + fn)
    # Calculate accuracy
    accuracy = (specificity * tnr + sensitivity) / 2
    return accuracy


def f1_at_specificity(y_true, y_pred, specificity, sample_weight=None):
    # Calculate the threshold corresponding to the specified specificity
    fpr, _, thresholds = roc_curve(
        y_true, y_pred, pos_label=1, sample_weight=sample_weight)
    idx = np.argmin(np.abs(fpr - (1 - specificity)))
    threshold = thresholds[idx]
    # Classify predictions based on the threshold
    y_pred_binary = y_pred > threshold
    # Calculate F1 score
    f1 = f1_score(y_true, y_pred_binary, sample_weight=sample_weight)
    return f1


model_metrics = [
    'accuracy',
    'auc'
]


def build_model(num_classes, model_input_shape):
    inputs = layers.Input(shape=model_input_shape)
    model = keras.applications.EfficientNetV2S(
        include_top=False,
        weights="imagenet",
        input_tensor=inputs,
        # classes=num_classes,
        # classifier_activation="sigmoid",
        include_preprocessing=True)
    # Freeze the pretrained weights
    model.trainable = False
    # Rebuild top
    x = layers.GlobalAveragePooling2D(name="avg_pool")(model.output)
    x = layers.BatchNormalization()(x)
    top_dropout_rate = 0.2
    x = layers.Dropout(top_dropout_rate, name="top_dropout")(x)
    outputs = layers.Dense(num_classes, activation="sigmoid", name="pred")(x)
    # Compile
    model = keras.Model(inputs, outputs, name="EfficientNet")
    optimizer = keras.optimizers.Adam(learning_rate=model_learning_rate)
    model.compile(
        optimizer=optimizer, loss="binary_crossentropy", metrics=model_metrics
    )
    return model


def unfreeze_model(model):
    # We unfreeze the top 20 layers while leaving BatchNorm layers frozen
    for layer in model.layers[-20:]:
        if not isinstance(layer, layers.BatchNormalization):
            layer.trainable = True
    optimizer = keras.optimizers.Adam(
        learning_rate=model_finetune_learning_rate)
    model.compile(
        optimizer=optimizer,
        loss="binary_crossentropy",
        metrics=model_metrics
    )


def split_data(x, y, z, fold, seed):
    kfold_grouped_for_validation = StratifiedGroupKFold(
        n_splits=fold, shuffle=True, random_state=seed)
    # get the first fold
    for i, (train_index, val_index) in enumerate(kfold_grouped_for_validation.split(x, y, groups=z)):
        x_train, x_val = x[train_index], x[val_index]
        y_train, y_val = y[train_index], y[val_index]
        if i == 0:
            break
    return x_train, x_val, y_train, y_val


def myprint(s):
    with open(roc_curve_csv_path_ind_test + model_name + '_model_architecture_summary.txt', 'a') as f:
        print(s, file=f)
        # report which file we are on
        print('model architecture summary saved to: ',
              roc_curve_csv_path_ind_test + model_name + '_model_architecture_summary.txt')


# if scale_pixel_to_255:
#     x = x * 255  # this is needed for efficientnet preprocessing
#     x_indtest = x_indtest * 255

if 'efficientnet' in model_name.lower() or 'efficient_net' in model_name.lower():
    x = x * 255  # this is needed for efficientnet preprocessing
    x_indtest = x_indtest * 255

# check if the shape of x and x_indtest are the same
if x.shape[1:] != x_indtest.shape[1:]:
    raise ValueError('x and x_indtest have different shapes')

# remove the problem bin if needed and x.shape[1] = 465
if remove_problem_bin and x.shape[1] == 465:
    x = np.delete(x, 276, axis=1)
    x_indtest = np.delete(x_indtest, 276, axis=1)

# build the model
model = build_model(num_classes, model_input_shape=x.shape[1:])
# save model summary to a file
model.summary(print_fn=myprint)

x_train, x_val, y_train, y_val = split_data(
    x, y, z, fold=train_val_split_fold, seed=seed)


model.fit(x_train,
          y_train,
          epochs=model_epochs,
          batch_size=model_batch_size,
          class_weight=model_class_weight,
          validation_data=(x_val, y_val))

# fine tune
unfreeze_model(model)
model.fit(x_train,
          y_train,
          epochs=model_fine_tune_epochs,
          batch_size=model_batch_size,
          class_weight=model_class_weight,
          validation_data=(x_val, y_val))

y_pred = model.predict(x_indtest)
y_pred_class = np.where(y_pred > 0.5, 1, 0)


###############################################################################
# metrics
###############################################################################
which_fold = 0
which_repeat = 0
tn, fp, fn, tp = confusion_matrix(y_indtest, y_pred_class).ravel()
print(confusion_matrix(y_indtest, y_pred_class))
print('tn: ', tn)
print('fp: ', fp)
print('fn: ', fn)
print('tp: ', tp)

# print classification report
print(classification_report(y_indtest, y_pred_class))

# make a dataframe including  cols, y_indtest, y_pred, z_test
y_pred_f = y_pred.flatten()
y_pred_class_f = y_pred_class.flatten()
y_test_f = y_indtest.flatten()
z_test_f = z_indtest.flatten()
bam_id_test_f = bam_id_indtest.flatten()

# create a dataframe to store y_test_f, y_pred_f, z_test_f
y_pred_df = pd.DataFrame(
    {'repeat': which_repeat,
     'fold': which_fold,
     'model_name': model_name,
     'y_test': y_test_f,
     'y_pred': y_pred_f,
     'y_pred_class': y_pred_class_f,
     'z_test': z_test_f,
     'bam_id_test': bam_id_test_f})

y_pred_df_filename = 'repeat_{}_fold_{}_y_test_y_pred_z_test.csv'.format(
    which_repeat, which_fold)
y_pred_df_filename_fullpath = roc_curve_csv_path_ind_test + y_pred_df_filename
y_pred_df.to_csv(y_pred_df_filename_fullpath, index=False)

print('y_pred_df file saved to: ', y_pred_df_filename_fullpath)

# if_set_sample_weight_for_testing = True, calculate the sample_weight for y_indtest
if if_set_sample_weight_for_testing:
    from sklearn.utils import class_weight
    sample_weight = class_weight.compute_sample_weight(
        class_weight='balanced', y=y_indtest)
else:
    sample_weight = None

# add repeat number
results['repeat'] = which_repeat
# add fold number
results['fold'] = which_fold

# add model_name
results['model_name'] = model_name
# calculate the the auroc
auroc = roc_auc_score(y_indtest, y_pred, sample_weight=sample_weight)
results['test_auroc'] = auroc
print('auroc: ', auroc)

# calculate the the auprc
auprc = average_precision_score(
    y_indtest, y_pred, sample_weight=sample_weight)
results['test_auprc'] = auprc
print('auprc: ', auprc)

# calculate generic accuracy
acc = accuracy_score(y_indtest, y_pred > 0.5, sample_weight=sample_weight)
results['test_acc'] = acc
print('acc: ', acc)

# calculate generic ppv
ppv = precision_score(
    y_indtest, y_pred > 0.5, sample_weight=sample_weight)
results['test_ppv'] = ppv
print('ppv: ', ppv)

# calculate generic sensitivity
sensitivity = recall_score(
    y_indtest, y_pred > 0.5, sample_weight=sample_weight)
results['test_sensitivity'] = sensitivity
print('sensitivity: ', sensitivity)

# calculate generic specificity
tn, fp, fn, tp = confusion_matrix(
    y_indtest, y_pred > 0.5, sample_weight=sample_weight).ravel()
specificity = tn / (tn + fp)
results['test_specificity'] = specificity
print('specificity: ', specificity)

# calculate the f1 score
f1 = f1_score(y_indtest, y_pred > 0.5, sample_weight=sample_weight)
results['test_f1'] = f1
print('f1: ', f1)

# get roc curve and store in csv, file name as repeat_fold_roc.csv
fpr, tpr, thresholds = roc_curve(
    y_indtest, y_pred, pos_label=1, sample_weight=sample_weight)
roc_df = pd.DataFrame(
    {'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds})
roc_filename = 'repeat_{}_fold_{}_roc.csv'.format(
    which_repeat, which_fold)
roc_filename_fullpath = roc_curve_csv_path_ind_test + roc_filename
roc_df.to_csv(roc_filename_fullpath, index=False)
print('roc curve saved to: ', roc_filename_fullpath)

print('* ' * 50)

# get the sensitivity at 99% specificity, and threshold
specificity = 0.99
sensitivity_99spe, threshold_99spe = sensitivity_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_sen_{}spe'.format(
        int(specificity * 100))] = sensitivity_99spe
results['test_threshold_{}spe'.format(
        int(specificity * 100))] = threshold_99spe
results['test_ppv_threshold_{}spe'.format(int(
        specificity * 100))] = ppv_at_specificity(y_indtest, y_pred, threshold_99spe, sample_weight)

acc_99spe = accuracy_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_acc_{}spe'.format(int(specificity * 100))] = acc_99spe

f1_99spe = f1_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_f1_{}spe'.format(int(specificity * 100))] = f1_99spe

print('sensitivity_99spe: ', sensitivity_99spe)
print('threshold_99spe: ', threshold_99spe)
print('ppv_threshold_99spe: ', ppv_at_specificity(
    y_indtest, y_pred, threshold_99spe, sample_weight))
print('acc_99spe: ', acc_99spe)
print('f1_99spe: ', f1_99spe)

# get the sample id of TP, TN, FP, FN predictions at 99% specificity
tn, fp, fn, tp = confusion_matrix(
    y_indtest, y_pred > threshold_99spe, sample_weight=sample_weight).ravel()

cm_99spe_filename = 'repeat_{}_fold_{}_cm_99spec.csv'.format(
    which_repeat, which_fold)
cm_99spe_filename_fullpath = roc_curve_csv_path_ind_test + cm_99spe_filename
cm_99spe_df = pd.DataFrame(
    {'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp}, index=[0])
cm_99spe_df.to_csv(
    cm_99spe_filename_fullpath, index=False)

# make TN_99spe, FP_99spe, FN_99spe, TP_99spe folders in roc_curve_csv_path_ind_test
os.makedirs(roc_curve_csv_path_ind_test + 'TN_99spe/', exist_ok=True)
os.makedirs(roc_curve_csv_path_ind_test + 'FP_99spe/', exist_ok=True)
os.makedirs(roc_curve_csv_path_ind_test + 'FN_99spe/', exist_ok=True)
os.makedirs(roc_curve_csv_path_ind_test + 'TP_99spe/', exist_ok=True)

# get sample id of TP, TN, FP, FN predictions at 99% specificity
tn_sample_id = np.where((y_test_f == 0) & (y_pred_f < threshold_99spe))
fp_sample_id = np.where((y_test_f == 0) & (y_pred_f > threshold_99spe))
fn_sample_id = np.where((y_test_f == 1) & (y_pred_f < threshold_99spe))
tp_sample_id = np.where((y_test_f == 1) & (y_pred_f > threshold_99spe))

# print the sample_id
print('tn_sample_id: ', tn_sample_id)
print('fp_sample_id: ', fp_sample_id)
print('fn_sample_id: ', fn_sample_id)
print('tp_sample_id: ', tp_sample_id)

tn_sample_namelist = bam_id_indtest[tn_sample_id]
fp_sample_namelist = bam_id_indtest[fp_sample_id]
fn_sample_namelist = bam_id_indtest[fn_sample_id]
tp_sample_namelist = bam_id_indtest[tp_sample_id]

# make a dataframe for each of the TP, TN, FP, FN predictions, colnames are TP, TN, FP, FN
tn_sample_df = pd.DataFrame(tn_sample_namelist, columns=['TN'])
fp_sample_df = pd.DataFrame(fp_sample_namelist, columns=['FP'])
fn_sample_df = pd.DataFrame(fn_sample_namelist, columns=['FN'])
tp_sample_df = pd.DataFrame(tp_sample_namelist, columns=['TP'])

# save the tn_sample_df, fp_sample_df, fn_sample_df, tp_sample_df to csv
tn_sample_filename = 'repeat_{}_fold_{}_TN_samples_99spec.csv'.format(
    which_repeat, which_fold)
fp_sample_filename = 'repeat_{}_fold_{}_FP_samples_99spec.csv'.format(
    which_repeat, which_fold)
fn_sample_filename = 'repeat_{}_fold_{}_FN_samples_99spec.csv'.format(
    which_repeat, which_fold)
tp_sample_filename = 'repeat_{}_fold_{}_TP_samples_99spec.csv'.format(
    which_repeat, which_fold)

tn_sample_filename_fullpath = roc_curve_csv_path_ind_test + \
    'TN_99spe/' + tn_sample_filename
fp_sample_filename_fullpath = roc_curve_csv_path_ind_test + \
    'FP_99spe/' + fp_sample_filename
fn_sample_filename_fullpath = roc_curve_csv_path_ind_test + \
    'FN_99spe/' + fn_sample_filename
tp_sample_filename_fullpath = roc_curve_csv_path_ind_test + \
    'TP_99spe/' + tp_sample_filename

tn_sample_df.to_csv(tn_sample_filename_fullpath, index=False)
fp_sample_df.to_csv(fp_sample_filename_fullpath, index=False)
fn_sample_df.to_csv(fn_sample_filename_fullpath, index=False)
tp_sample_df.to_csv(tp_sample_filename_fullpath, index=False)
# --------------------------------------------------------------------------------
# get the sensitivity at 98% specificity, and threshold
# --------------------------------------------------------------------------------
specificity = 0.98
sensitivity_98spe, threshold_98spe = sensitivity_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_sen_{}spe'.format(
        int(specificity * 100))] = sensitivity_98spe
results['test_threshold_{}spe'.format(
        int(specificity * 100))] = threshold_98spe
results['test_ppv_threshold_{}spe'.format(int(
        specificity * 100))] = ppv_at_specificity(y_indtest, y_pred, threshold_98spe, sample_weight)
acc_98spe = accuracy_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_acc_{}spe'.format(int(specificity * 100))] = acc_98spe
f1_98spe = f1_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_f1_{}spe'.format(int(specificity * 100))] = f1_98spe

print('sensitivity_98spe: ', sensitivity_98spe)
print('threshold_98spe: ', threshold_98spe)
print('ppv_threshold_98spe: ', ppv_at_specificity(
    y_indtest, y_pred, threshold_98spe, sample_weight))
print('acc_98spe: ', acc_98spe)
print('f1_98spe: ', f1_98spe)

# -------------------------------------------------------------------------------
# get the sensitivity at 95% specificity, and threshold
# -------------------------------------------------------------------------------
specificity = 0.95
sensitivity_95spe, threshold_95spe = sensitivity_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_sen_{}spe'.format(
        int(specificity * 100))] = sensitivity_95spe
results['test_threshold_{}spe'.format(
        int(specificity * 100))] = threshold_95spe
results['test_ppv_threshold_{}spe'.format(int(
        specificity * 100))] = ppv_at_specificity(y_indtest, y_pred, threshold_95spe, sample_weight)
acc_95spe = accuracy_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_acc_{}spe'.format(int(specificity * 100))] = acc_95spe
f1_95spe = f1_at_specificity(
    y_indtest, y_pred, specificity, sample_weight)
results['test_f1_{}spe'.format(int(specificity * 100))] = f1_95spe

print('sensitivity_95spe: ', sensitivity_95spe)
print('threshold_95spe: ', threshold_95spe)
print('ppv_threshold_95spe: ', ppv_at_specificity(
    y_indtest, y_pred, threshold_95spe, sample_weight))
print('acc_95spe: ', acc_95spe)
print('f1_95spe: ', f1_95spe)

# -------------------------------------------------------------------------------
# Append values to the lists associated with each key
# -------------------------------------------------------------------------------
for key, value in results.items():
    results_container[key].append(value)
results_dict = dict(results_container)
results_df = pd.DataFrame.from_dict(results_dict)
# convet to dataframe
metrics_filename = 'repeat_{}_fold_{}_metrics.csv'.format(
    which_repeat, which_fold)
metrics_filename_fullpath = roc_curve_csv_path_ind_test + metrics_filename
results_df.to_csv(metrics_filename_fullpath, index=False)
print('metrics saved to: ', metrics_filename_fullpath)
print(results_df)
