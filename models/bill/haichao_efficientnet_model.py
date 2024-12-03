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

# parameters

parser = argparse.ArgumentParser(description='Script parameters')
parser.add_argument('--input_data_path', type=str,
                    default='/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn7/plasma_tp1_0.1x_model_C2/', help='Path to input data')
parser.add_argument('--tf_label', type=str,
                    default='0.03_0.1', help='TF label')
parser.add_argument('--model_name', type=str,
                    default='efficient_net_s', help='Model name')
parser.add_argument('--remove_problem_bin',
                    action='store_true', help='Remove problem bin')
parser.add_argument('--output_path', type=str,
                    default='/home/nrlab/wang04/ulyses/models/bill/', help='Path to output files')

parser.add_argument('--scale_pixel_to_255',
                    action='store_true', help='Scale pixel values to 255')
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
output_path = args.output_path
scale_pixel_to_255 = args.scale_pixel_to_255
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


# Set the random seed for reproducibility
np.random.seed(seed)
tf.random.set_seed(seed)
random.seed(seed)

# x_file_path = input_data_path + \
#     'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_{}_.npy'.format(tf_label)
# y_file_path = input_data_path + \
#     'y_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_{}_.npy'.format(tf_label)
# z_file_path = input_data_path + \
#     'z_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_{}_.npy'.format(tf_label)


# list all files in the input_data_path, find the x_file as file start with 'x_clean', and contains 'TF_{}'.format(tf_label), and ends with '.npy'
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

# create the roc_curve_csv_path if not exist
os.makedirs(roc_curve_csv_path, exist_ok=True)


# load the data
x = np.load(x_file_path, allow_pickle=True)
y = np.load(y_file_path, allow_pickle=True)
z = np.load(z_file_path, allow_pickle=True)
bam_id = np.load(bam_id_file_path, allow_pickle=True)

if scale_pixel_to_255:
    x = x * 255  # this is needed for efficientnet preprocessing
if remove_problem_bin and x.shape[1] == 465:
    x = np.delete(x, 276, axis=1)


model_input_shape = x.shape[1:]


# print the parameters and file paths
print('= ' * 50)
print('input_data_path: ', input_data_path)
print('tf_label: ', tf_label)
print('model_name: ', model_name)
print('remove_problem_bin: ', remove_problem_bin)
print('output_path: ', output_path)
print('scale_pixel_to_255: ', scale_pixel_to_255)
print('num_classes: ', num_classes)
print('model_learning_rate: ', model_learning_rate)
print('model_finetune_learning_rate: ', model_finetune_learning_rate)
print('num_repeats: ', num_repeats)
print('num_folds: ', num_folds)
print('train_val_split_fold: ', train_val_split_fold)
print('model_epochs: ', model_epochs)
print('model_batch_size: ', model_batch_size)
print('model_fine_tune_epochs: ', model_fine_tune_epochs)
print('set_sample_weight_for_testing: ', if_set_sample_weight_for_testing)
print('x_file_path: ', x_file_path)
print('y_file_path: ', y_file_path)
print('z_file_path: ', z_file_path)
print('performance_csv_path: ', performance_csv_path)
print('roc_curve_csv_path: ', roc_curve_csv_path)
print('seed', seed)
print('x shape: ', x.shape)
print('y shape: ', y.shape)
print('z shape: ', z.shape)
print('bam_id shape: ', bam_id.shape)


# save parameters to a txt file
with open(parameters_csv_path, 'w') as f:
    print('input_data_path: ', input_data_path, file=f)
    print('tf_label: ', tf_label, file=f)
    print('model_name: ', model_name, file=f)
    print('remove_problem_bin: ', remove_problem_bin, file=f)
    print('output_path: ', output_path, file=f)
    print('scale_pixel_to_255: ', scale_pixel_to_255, file=f)
    print('num_classes: ', num_classes, file=f)
    print('model_learning_rate: ', model_learning_rate, file=f)
    print('model_finetune_learning_rate: ',
          model_finetune_learning_rate, file=f)
    print('num_repeats: ', num_repeats, file=f)
    print('num_folds: ', num_folds, file=f)
    print('train_val_split_fold: ', train_val_split_fold, file=f)
    print('model_epochs: ', model_epochs, file=f)
    print('model_batch_size: ', model_batch_size, file=f)
    print('model_fine_tune_epochs: ', model_fine_tune_epochs, file=f)
    print('set_sample_weight_for_testing: ',
          if_set_sample_weight_for_testing, file=f)
    print('x_file_path: ', x_file_path, file=f)
    print('y_file_path: ', y_file_path, file=f)
    print('z_file_path: ', z_file_path, file=f)
    print('performance_csv_path: ', performance_csv_path, file=f)
    print('roc_curve_csv_path: ', roc_curve_csv_path, file=f)
    print('seed', seed, file=f)
    print('x shape: ', x.shape, file=f)
    print('y shape: ', y.shape, file=f)
    print('z shape: ', z.shape, file=f)
    print('bam_id shape: ', bam_id.shape, file=f)
    print('parameters settings saved to: ',
          parameters_csv_path)

print('parameters saved to: ', parameters_csv_path)
print('= ' * 50)
# prepare the data


# compute class weights
model_class_weight = dict(
    zip(np.unique(y),
        class_weight.compute_class_weight(class_weight='balanced',
                                          classes=np.unique(np.array(y)), y=y)))


# def get_f1(y_true, y_pred):  # taken from old keras source code
#    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
#    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
#    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
#    precision = true_positives / (predicted_positives + K.epsilon())
#    recall = true_positives / (possible_positives + K.epsilon())
#    f1_val = 2*(precision*recall)/(precision+recall+K.epsilon())
#    return f1_val


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
    tf.keras.metrics.AUC(curve="ROC", name='AUROC'),
    tf.keras.metrics.AUC(curve="PR", name='AUPRC'),
    tf.keras.metrics.Precision(name='PPV'),
    tf.keras.metrics.Recall(name='recall_sens'),
    # add f1 score
    # get_f1,
    # keras.metrics.F1Score(name='f1'),
    tf.keras.metrics.SensitivityAtSpecificity(
        specificity=0.95, name="sen_95spe"),
    tf.keras.metrics.SensitivityAtSpecificity(
        specificity=0.98, name="sen_98spe"),
    tf.keras.metrics.SensitivityAtSpecificity(
        specificity=0.99, name="sen_99spe")
]

model_metrics = [
    'accuracy',
    'auc'
]


def build_model(num_classes):
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


# build the model
model_tmp = build_model(num_classes)


def myprint(s):
    with open(output_path + model_name + '_model_architecture_summary.txt', 'a') as f:
        print(s, file=f)
        # report which file we are on
        print('model architecture summary saved to: ',
              output_path + model_name + '_model_architecture_summary.txt')


model_tmp.summary(print_fn=myprint)


###############################################################################
# REPEAT
###############################################################################
results = {}
results_container = defaultdict(list)

# set random seed
np.random.seed(seed)

for i in range(num_repeats):
    # Generate a new seed for each repeat
    repeat_seed = seed + i
    np.random.seed(repeat_seed)
    tf.random.set_seed(repeat_seed)
    random.seed(repeat_seed)

    # report which loop we are on
    print('= ' * 50)
    print('Repeat: ', i + 1)
    print('Seed: ', repeat_seed)

    kfold_grouped = StratifiedGroupKFold(
        n_splits=num_folds, shuffle=True, random_state=repeat_seed)

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

    # split into train and test
    for fold_i, (train_idx, test_idx) in enumerate(kfold_grouped.split(x, y, groups=z)):
        # print repeat and fold number
        print('- ' * 50)
        print('Repeat: ', i + 1)
        print('Fold: ', fold_i + 1)
        print('Seed: ', repeat_seed)

        # Split the data into train and test sets
        x_train, x_test = x[train_idx], x[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        z_train, z_test = z[train_idx], z[test_idx]
        bam_id_train, bam_id_test = bam_id[train_idx], bam_id[test_idx]

        # manually select validation data using grouped stratified way
        x_train, x_val, y_train, y_val = split_data(
            x_train, y_train, z_train, fold=train_val_split_fold, seed=repeat_seed)

        # Create and train the model
        model = build_model(num_classes)
        # print model
        # model.summary()

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

        # using test data
        print("Evaluate...")
        # Store the performance metrics for each repetition

        # model prediction
        y_pred = model.predict(x_test)
        y_pred_class = np.where(y_pred > 0.5, 1, 0)

        # # use grad-CAM to show the heatmap of one sample
        # # get the last convolutional layer
        # last_conv_layer = model.get_layer('top_activation')
        # # get the gradient tape
        # grad_model = tf.keras.models.Model(
        #     [model.inputs], [last_conv_layer.output, model.output])

        # make a dataframe including 3 cols, y_test, y_pred, z_test
        y_pred_f = y_pred.flatten()
        y_pred_class_f = y_pred_class.flatten()
        y_test_f = y_test.flatten()
        z_test_f = z_test.flatten()
        bam_id_test_f = bam_id_test.flatten()

        # create a dataframe to store y_test_f, y_pred_f, z_test_f
        y_pred_df = pd.DataFrame(
            {'repeat': i + 1,
             'fold': fold_i + 1,
             'model_name': model_name,
             'y_test': y_test_f,
             'y_pred': y_pred_f,
             'y_pred_class': y_pred_class_f,
             'z_test': z_test_f,
             'bam_id_test': bam_id_test_f})

        y_pred_df_filename = 'repeat_{}_fold_{}_y_test_y_pred_z_test.csv'.format(
            i + 1, fold_i + 1)
        y_pred_df_filename_fullpath = roc_curve_csv_path + y_pred_df_filename
        y_pred_df.to_csv(y_pred_df_filename_fullpath, index=False)

        # if_set_sample_weight_for_testing = True, calculate the sample_weight for y_test
        if if_set_sample_weight_for_testing:
            from sklearn.utils import class_weight
            sample_weight = class_weight.compute_sample_weight(
                class_weight='balanced', y=y_test)
        else:
            sample_weight = None

        # add repeat number
        results['repeat'] = i + 1
        # add fold number
        results['fold'] = fold_i + 1

        # add model name
        results['model_name'] = model_name

        # calculate the the auroc
        from sklearn.metrics import roc_auc_score
        auroc = roc_auc_score(y_test, y_pred, sample_weight=sample_weight)
        results['test_auroc'] = auroc
        print('auroc: ', auroc)

        # calculate the the auprc
        from sklearn.metrics import average_precision_score
        auprc = average_precision_score(
            y_test, y_pred, sample_weight=sample_weight)
        results['test_auprc'] = auprc
        print('auprc: ', auprc)

        # calculate generic accuracy
        from sklearn.metrics import accuracy_score
        acc = accuracy_score(y_test, y_pred > 0.5, sample_weight=sample_weight)
        results['test_acc'] = acc
        print('acc: ', acc)

        # calculate generic ppv
        ppv = precision_score(
            y_test, y_pred > 0.5, sample_weight=sample_weight)
        results['test_ppv'] = ppv
        print('ppv: ', ppv)

        # calculate generic sensitivity
        from sklearn.metrics import recall_score
        sensitivity = recall_score(
            y_test, y_pred > 0.5, sample_weight=sample_weight)
        results['test_sensitivity'] = sensitivity
        print('sensitivity: ', sensitivity)

        # calculate generic specificity
        from sklearn.metrics import confusion_matrix
        tn, fp, fn, tp = confusion_matrix(
            y_test, y_pred > 0.5, sample_weight=sample_weight).ravel()
        specificity = tn / (tn + fp)
        results['test_specificity'] = specificity
        print('specificity: ', specificity)

        # calculate the f1 score
        f1 = f1_score(y_test, y_pred > 0.5, sample_weight=sample_weight)
        results['test_f1'] = f1
        print('f1: ', f1)

        # get roc curve and store in csv, file name as repeat_fold_roc.csv
        fpr, tpr, thresholds = roc_curve(
            y_test, y_pred, pos_label=1, sample_weight=sample_weight)
        roc_df = pd.DataFrame(
            {'fpr': fpr, 'tpr': tpr, 'thresholds': thresholds})
        roc_filename = 'repeat_{}_fold_{}_roc.csv'.format(i + 1, fold_i + 1)
        roc_filename_fullpath = roc_curve_csv_path + roc_filename
        roc_df.to_csv(roc_filename_fullpath, index=False)
        print('roc curve saved to: ', roc_filename_fullpath)

        print('* ' * 50)

        # get the sensitivity at 99% specificity, and threshold
        specificity = 0.99
        sensitivity_99spe, threshold_99spe = sensitivity_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_sen_{}spe'.format(
            int(specificity * 100))] = sensitivity_99spe
        results['test_threshold_{}spe'.format(
            int(specificity * 100))] = threshold_99spe
        results['test_ppv_threshold_{}spe'.format(int(
            specificity * 100))] = ppv_at_specificity(y_test, y_pred, threshold_99spe, sample_weight)

        acc_99spe = accuracy_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_acc_{}spe'.format(int(specificity * 100))] = acc_99spe

        f1_99spe = f1_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_f1_{}spe'.format(int(specificity * 100))] = f1_99spe

        print('sensitivity_99spe: ', sensitivity_99spe)
        print('threshold_99spe: ', threshold_99spe)
        print('ppv_threshold_99spe: ', ppv_at_specificity(
            y_test, y_pred, threshold_99spe, sample_weight))
        print('acc_99spe: ', acc_99spe)
        print('f1_99spe: ', f1_99spe)

        # get the sample id of TP, TN, FP, FN predictions at 99% specificity
        tn, fp, fn, tp = confusion_matrix(
            y_test, y_pred > threshold_99spe, sample_weight=sample_weight).ravel()

        cm_99spe_filename = 'repeat_{}_fold_{}_cm_99spec.csv'.format(
            i + 1, fold_i + 1)
        cm_99spe_filename_fullpath = roc_curve_csv_path + cm_99spe_filename
        cm_99spe_df = pd.DataFrame(
            {'tn': tn, 'fp': fp, 'fn': fn, 'tp': tp}, index=[0])
        cm_99spe_df.to_csv(
            cm_99spe_filename_fullpath, index=False)

        # make TN_99spe, FP_99spe, FN_99spe, TP_99spe folders in roc_curve_csv_path
        os.makedirs(roc_curve_csv_path + 'TN_99spe/', exist_ok=True)
        os.makedirs(roc_curve_csv_path + 'FP_99spe/', exist_ok=True)
        os.makedirs(roc_curve_csv_path + 'FN_99spe/', exist_ok=True)
        os.makedirs(roc_curve_csv_path + 'TP_99spe/', exist_ok=True)

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

        tn_sample_namelist = z_test[tn_sample_id]
        fp_sample_namelist = z_test[fp_sample_id]
        fn_sample_namelist = z_test[fn_sample_id]
        tp_sample_namelist = z_test[tp_sample_id]

        # make a dataframe for each of the TP, TN, FP, FN predictions, colnames are TP, TN, FP, FN
        tn_sample_df = pd.DataFrame(tn_sample_namelist, columns=['TN'])
        fp_sample_df = pd.DataFrame(fp_sample_namelist, columns=['FP'])
        fn_sample_df = pd.DataFrame(fn_sample_namelist, columns=['FN'])
        tp_sample_df = pd.DataFrame(tp_sample_namelist, columns=['TP'])

        # save the tn_sample_df, fp_sample_df, fn_sample_df, tp_sample_df to csv
        tn_sample_filename = 'repeat_{}_fold_{}_TN_samples_99spec.csv'.format(
            i + 1, fold_i + 1)
        fp_sample_filename = 'repeat_{}_fold_{}_FP_samples_99spec.csv'.format(
            i + 1, fold_i + 1)
        fn_sample_filename = 'repeat_{}_fold_{}_FN_samples_99spec.csv'.format(
            i + 1, fold_i + 1)
        tp_sample_filename = 'repeat_{}_fold_{}_TP_samples_99spec.csv'.format(
            i + 1, fold_i + 1)

        tn_sample_filename_fullpath = roc_curve_csv_path + 'TN_99spe/' + tn_sample_filename
        fp_sample_filename_fullpath = roc_curve_csv_path + 'FP_99spe/' + fp_sample_filename
        fn_sample_filename_fullpath = roc_curve_csv_path + 'FN_99spe/' + fn_sample_filename
        tp_sample_filename_fullpath = roc_curve_csv_path + 'TP_99spe/' + tp_sample_filename

        tn_sample_df.to_csv(tn_sample_filename_fullpath, index=False)
        fp_sample_df.to_csv(fp_sample_filename_fullpath, index=False)
        fn_sample_df.to_csv(fn_sample_filename_fullpath, index=False)
        tp_sample_df.to_csv(tp_sample_filename_fullpath, index=False)

        # get the sensitivity at 98% specificity, and threshold
        specificity = 0.98
        sensitivity_98spe, threshold_98spe = sensitivity_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_sen_{}spe'.format(
            int(specificity * 100))] = sensitivity_98spe
        results['test_threshold_{}spe'.format(
            int(specificity * 100))] = threshold_98spe
        results['test_ppv_threshold_{}spe'.format(int(
            specificity * 100))] = ppv_at_specificity(y_test, y_pred, threshold_98spe, sample_weight)
        acc_98spe = accuracy_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_acc_{}spe'.format(int(specificity * 100))] = acc_98spe
        f1_98spe = f1_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_f1_{}spe'.format(int(specificity * 100))] = f1_98spe

        print('sensitivity_98spe: ', sensitivity_98spe)
        print('threshold_98spe: ', threshold_98spe)
        print('ppv_threshold_98spe: ', ppv_at_specificity(
            y_test, y_pred, threshold_98spe, sample_weight))
        print('acc_98spe: ', acc_98spe)
        print('f1_98spe: ', f1_98spe)

        # get the sensitivity at 95% specificity, and threshold
        specificity = 0.95
        sensitivity_95spe, threshold_95spe = sensitivity_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_sen_{}spe'.format(
            int(specificity * 100))] = sensitivity_95spe
        results['test_threshold_{}spe'.format(
            int(specificity * 100))] = threshold_95spe
        results['test_ppv_threshold_{}spe'.format(int(
            specificity * 100))] = ppv_at_specificity(y_test, y_pred, threshold_95spe, sample_weight)
        acc_95spe = accuracy_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_acc_{}spe'.format(int(specificity * 100))] = acc_95spe
        f1_95spe = f1_at_specificity(
            y_test, y_pred, specificity, sample_weight)
        results['test_f1_{}spe'.format(int(specificity * 100))] = f1_95spe

        print('sensitivity_95spe: ', sensitivity_95spe)
        print('threshold_95spe: ', threshold_95spe)
        print('ppv_threshold_95spe: ', ppv_at_specificity(
            y_test, y_pred, threshold_95spe, sample_weight))
        print('acc_95spe: ', acc_95spe)
        print('f1_95spe: ', f1_95spe)

        # Append values to the lists associated with each key
        for key, value in results.items():
            results_container[key].append(value)
        # convert results to a dictionary
        results_dict = dict(results_container)
        results_df = pd.DataFrame.from_dict(results_dict)
        # convet to dataframe
        metrics_filename = 'repeat_{}_fold_{}_metrics.csv'.format(
            i + 1, fold_i + 1)
        metrics_filename_fullpath = roc_curve_csv_path + metrics_filename
        results_df.to_csv(metrics_filename_fullpath, index=False)
        print('metrics saved to: ', metrics_filename_fullpath)

        print(results_df)


###############################################################################
# store results in a csv
###############################################################################

# Create a pandas DataFrame from the results dictionary
# remove 'estimator' key
relevant_keys = [k for k in results_container.keys() if k != 'estimator']
filtered_dict = {k: v for k, v in results_container.items()
                 if k in relevant_keys}

# create a pandas DataFrame from the new dictionary
cv_df = pd.DataFrame.from_dict(filtered_dict)
# save to csv
print(cv_df)
cv_df.to_csv(performance_csv_path, index=False)
print('Results saved to: ', performance_csv_path)

# clean up memory
del results
tf.keras.backend.clear_session()
