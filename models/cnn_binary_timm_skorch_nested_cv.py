import pickle
from sklearn.metrics import classification_report
from prettytable import PrettyTable
from skorch.callbacks.scoring import EpochScoring
from skorch.callbacks.logging import SacredLogger
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score
import numpy as np
import pandas as pd
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
# grid search and random search
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
# torch
import torch
import torch.nn as nn
import torch.optim as optim
from torch import nn
import timm
import fastai

# skorch
from skorch.callbacks import ProgressBar
from skorch.callbacks import Freezer
from skorch.callbacks import Checkpoint
from skorch.callbacks import LRScheduler
from skorch import NeuralNetClassifier
import skorch
from skorch import NeuralNetClassifier, NeuralNetBinaryClassifier
from skorch.helper import predefined_split
from skorch.dataset import Dataset


# parameters

parser = argparse.ArgumentParser(description='Script parameters')
parser.add_argument('--input_data_path', type=str,
                    default='/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/cnn_model/cnn11/plasma_tp1_0.1x_model_C5', help='Path to input data')
parser.add_argument('--tf_label', type=str,
                    default='0_0.03', help='TF label')
parser.add_argument('--model_name', type=str,
                    default='resnet34', help='Model name')
parser.add_argument('--remove_problem_bin',
                    action='store_true', help='Remove problem bin')
parser.add_argument('--num_classes', type=int,
                    default=1, help='Number of classes')
parser.add_argument('--model_learning_rate', type=float,
                    default=0.005, help='Learning rate for the model')
parser.add_argument('--model_finetune_learning_rate', type=float,
                    default=1e-5, help='Learning rate for fine-tuning the model')
parser.add_argument('--num_repeats', type=int,
                    default=10, help='Number of repeats')
parser.add_argument('--outer_cv_fold', type=int, default=5,
                    help='Number of folds for cross-validation')
parser.add_argument('--train_val_split_fold', type=int,
                    default=3, help='Fold number for train-validation split')
parser.add_argument('--model_epochs', type=int, default=30,
                    help='Number of epochs for training the model')
parser.add_argument('--model_batch_size', type=int, default=30,
                    help='Batch size for training the model')
parser.add_argument('--model_fine_tune_epochs', type=int,
                    default=5, help='Number of epochs for fine-tuning the model')
parser.add_argument('--set_sample_weight_for_testing',
                    action='store_true', help='Set sample weight for testing')
# train_from_scratch
parser.add_argument('--train_from_scratch',
                    action='store_true', help='Train the model from scratch')

# add a parameter called "seed" to the parser
parser.add_argument('--seed', type=int, default=0, help='Random seed')

# add a param called "hyperparameter_search_strategy", with default to "grid", choices=["grid", "random"]
parser.add_argument('--hyperparameter_search_strategy', type=str, default='grid',
                    choices=['grid', 'random'], help='Hyperparameter search strategy')
# add a param called "random_search_iter", with default to 100
parser.add_argument('--random_search_iter', type=int,
                    default=10, help='Number of iterations for random search')
parser.add_argument('--inner_cv_fold', type=int, default=3,
                    help='Number of inner cross-validation folds')
# add a param called "top_layers", default to "new", choices=["new", "original"]
parser.add_argument('--top_layers', type=str, default='new',
                    choices=['new', 'original'], help='Top layers of the model')

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
outer_cv_fold = args.outer_cv_fold
train_val_split_fold = args.train_val_split_fold
model_epochs = args.model_epochs
model_batch_size = args.model_batch_size
model_fine_tune_epochs = args.model_fine_tune_epochs
if_set_sample_weight_for_testing = args.set_sample_weight_for_testing
seed = args.seed
train_from_scratch = args.train_from_scratch
hyperparameter_search_strategy = args.hyperparameter_search_strategy
random_search_iter = args.random_search_iter
inner_cv_fold = args.inner_cv_fold
top_layers = args.top_layers


output_path = input_data_path
# Set the random seed for reproducibility
np.random.seed(seed)
random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

x_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'x_clean_layers') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
y_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'y_clean_layers') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
z_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'z_clean_layers') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]
bam_id_file_path_file = [f for f in os.listdir(input_data_path) if f.startswith(
    'bam_id_clean_layers') and 'TF_{}'.format(tf_label) in f and f.endswith('.npy')][0]

x_file_path = os.path.join(input_data_path, x_file_path_file)
y_file_path = os.path.join(input_data_path, y_file_path_file)
z_file_path = os.path.join(input_data_path, z_file_path_file)
bam_id_file_path = os.path.join(input_data_path, bam_id_file_path_file)

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
roc_curve_csv_path_ind_test = roc_curve_csv_path + \
    'nested_cv_' + hyperparameter_search_strategy + '/' + model_name + '/'

# create the roc_curve_csv_path if not exist
os.makedirs(roc_curve_csv_path, exist_ok=True)
os.makedirs(roc_curve_csv_path_ind_test, exist_ok=True)

# load the data
x = np.load(x_file_path, allow_pickle=True)
y = np.load(y_file_path, allow_pickle=True)
z = np.load(z_file_path, allow_pickle=True)
bam_id = np.load(bam_id_file_path, allow_pickle=True)

# load independent test data
x_test = np.load(x_indtest_file_path, allow_pickle=True)
y_test = np.load(y_indtest_file_path, allow_pickle=True)
z_test = np.load(z_indtest_file_path, allow_pickle=True)
bam_id_test = np.load(bam_id_indtest_file_path, allow_pickle=True)

# stop if the first dimension of x, y, z, bam_id are not the same
if x.shape[0] != y.shape[0] or x.shape[0] != z.shape[0] or x.shape[0] != bam_id.shape[0]:
    raise ValueError(
        'The first dimension of x, y, z, bam_id are not the same')

# stop if the first dimension of x_test, y_test, z_test, bam_id_test are not the same
if x_test.shape[0] != y_test.shape[0] or x_test.shape[0] != z_test.shape[0] or x_test.shape[0] != bam_id_test.shape[0]:
    raise ValueError(
        'The first dimension of x_test, y_test, z_test, bam_id_test are not the same')


# for x and x_test, if there are only 3 dimensions, add a new dimension as forth dimension
if len(x.shape) == 3:
    x = x[:, :, :, np.newaxis]
if len(x_test.shape) == 3:
    x_test = x_test[:, :, :, np.newaxis]


# if scale_pixel_to_255:
#     x = x * 255  # this is needed for efficientnet preprocessing
#     x_test = x_test * 255
# remove the problem bin if needed and x.shape[1] = 465
if remove_problem_bin and x.shape[1] == 465:
    x = np.delete(x, 276, axis=1)
    x_test = np.delete(x_test, 276, axis=1)

model_input_shape = x.shape[1:]

# get the number of input channels, which is the last dimension of the model_input_shape
n_input_channel = model_input_shape[-1]


# print the parameters and file paths
print('= ' * 50)
print('input_data_path: ', input_data_path)
print('tf_label: ', tf_label)
print('model_name: ', model_name)
print('remove_problem_bin: ', remove_problem_bin)
print('output_path: ', output_path)
# print('scale_pixel_to_255: ', scale_pixel_to_255)
print('num_classes: ', num_classes)
print('model_learning_rate: ', model_learning_rate)
print('model_finetune_learning_rate: ', model_finetune_learning_rate)
print('num_repeats: ', num_repeats)
print('outer_cv_fold: ', outer_cv_fold)
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
print("hyperparameter_search_strategy: ", hyperparameter_search_strategy)
print("random_search_iter: ", random_search_iter)
print("inner_cv_fold: ", inner_cv_fold)
print("train_from_scratch: ", train_from_scratch)
print("top_layers: ", top_layers)

print('= ' * 50)


# save parameters to a txt file
with open(parameters_csv_path, 'w') as f:
    print('input_data_path: ', input_data_path, file=f)
    print('tf_label: ', tf_label, file=f)
    print('model_name: ', model_name, file=f)
    print('remove_problem_bin: ', remove_problem_bin, file=f)
    print('output_path: ', output_path, file=f)
    # print('scale_pixel_to_255: ', scale_pixel_to_255, file=f)
    print('num_classes: ', num_classes, file=f)
    print('model_learning_rate: ', model_learning_rate, file=f)
    print('model_finetune_learning_rate: ',
          model_finetune_learning_rate, file=f)
    print('num_repeats: ', num_repeats, file=f)
    print('outer_cv_fold: ', outer_cv_fold, file=f)
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
    print("hyperparameter_search_strategy: ",
          hyperparameter_search_strategy, file=f)
    print("random_search_iter: ", random_search_iter, file=f)
    print("inner_cv_fold: ", inner_cv_fold, file=f)
    print("train_from_scratch: ", train_from_scratch, file=f)
    print("top_layers: ", top_layers, file=f)


print('parameters saved to: ', parameters_csv_path)
print('= ' * 50)
# prepare the data


# compute class weights
model_class_weight = dict(
    zip(np.unique(y),
        class_weight.compute_class_weight(class_weight='balanced',
                                          classes=np.unique(np.array(y)), y=y)))


# create a StratifiedGroupKFold object
def sensitivity_at_specificity(y_true, y_pred, specificity, sample_weight=None):
    fpr, tpr, thresholds = roc_curve(
        y_true, y_pred, pos_label=1, sample_weight=sample_weight)
    idx = np.argmin(np.abs(fpr - (1 - specificity)))
    return tpr[idx], thresholds[idx]


def ppv_at_specificity(y_true, y_pred, threshold, sample_weight=None):
    y_pred_above_threshold = y_pred > threshold
    y_pred_above_threshold = y_pred_above_threshold.astype(int)
    return precision_score(y_true, y_pred_above_threshold,
                           sample_weight=sample_weight,
                           zero_division=0)

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


def count_trainable_parameters(model):
    table = PrettyTable(["Modules", "Parameters"])
    total_params = 0
    for name, parameter in model.named_parameters():
        if not parameter.requires_grad:
            continue
        params = parameter.numel()
        table.add_row([name, params])
        total_params += params
    print(table)
    print(f"Total Trainable Params: {total_params}")
    return total_params


def count_total_parameters(model):
    table = PrettyTable(["Modules", "Parameters"])
    total_params = 0
    for name, parameter in model.named_parameters():
        params = parameter.numel()
        table.add_row([name, params])
        total_params += params
    print(table)
    print(f"Total Params: {total_params}")
    return total_params


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


results = {}
results_container = defaultdict(list)
# set the random seed for reproducibility

# scale the pixel values to 255 if model_name contains "efficientnet"
if 'efficientnet' in model_name.lower() or 'efficient_net' in model_name.lower():
    x = x * 255
    x_test = x_test * 255


# where the loop enters
for i in range(num_repeats):
    repeat_seed = seed + i
    np.random.seed(repeat_seed)
    random.seed(repeat_seed)
    which_repeat = i + 1
    print('= ' * 50)
    print('Repeat: ', which_repeat)
    print('Seed: ', repeat_seed)

    outer_cv = StratifiedGroupKFold(
        n_splits=outer_cv_fold, shuffle=True, random_state=repeat_seed)
    # randomized grid search
    inner_cv = StratifiedGroupKFold(
        n_splits=inner_cv_fold, shuffle=True, random_state=repeat_seed)

    model = timm.create_model(model_name=model_name,
                              pretrained=True,
                              in_chans=n_input_channel,
                              num_classes=1)

    # split into train and test
    for fold_i, (train_idx, test_idx) in enumerate(outer_cv.split(x, y, groups=z)):
        # print repeat and fold number
        which_fold = fold_i + 1
        print('- ' * 50)
        print('Repeat: ', which_repeat)
        print('Fold: ', which_fold)
        print('Seed: ', repeat_seed)
        # Split the data into train and test sets
        x_train, x_test = x[train_idx], x[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        z_train, z_test = z[train_idx], z[test_idx]
        bam_id_train, bam_id_test = bam_id[train_idx], bam_id[test_idx]
        within_fold_classes_weights = class_weight.compute_sample_weight(
            class_weight='balanced',
            y=y_train
        )

        if top_layers == 'new':
            num_in_features = model.get_classifier().in_features
            if 'efficientnet' in model_name.lower() or 'efficient_net' in model_name.lower():
                model.classifier = nn.Sequential(
                    # add batch norm
                    nn.BatchNorm1d(num_in_features),
                    nn.Dropout(0.2),
                    nn.Linear(num_in_features, out_features=1, bias=True)
                )
            elif 'resnet' in model_name.lower():
                model.fc = nn.Sequential(
                    # add batch norm
                    nn.BatchNorm1d(num_in_features),
                    nn.Dropout(0.2),
                    nn.Linear(num_in_features, out_features=1, bias=True)
                )
            else:
                model.classifier = nn.Sequential(
                    # add batch norm
                    nn.BatchNorm1d(num_in_features),
                    nn.Dropout(0.2),
                    nn.Linear(num_in_features, out_features=1, bias=True)
                )

        params = {
            # need to add more lower learning rates,e.g., 0.0001, 0.0003, 0.0005
            "lr": [0.0001, 0.0005, 0.001, 0.003, 0.005],
            "max_epochs": [10, 30, 40],
            "batch_size": [30, 40, 50, 60, 80]
        }

        train_class_weight = torch.Tensor(
            class_weight.compute_class_weight(class_weight='balanced', classes=np.unique(y), y=y))

        print('sample_weight: ', train_class_weight)

        pos_weight = torch.Tensor([train_class_weight[1] / train_class_weight[0]])
        print('pos_weight: ', pos_weight)

        # count_parameters(model)
        if train_from_scratch:
            net = NeuralNetBinaryClassifier(
                model,
                max_epochs=model_epochs,
                lr=model_learning_rate,
                batch_size=model_batch_size,
                criterion=nn.BCEWithLogitsLoss(),
                optimizer=optim.Adam,
                device='cpu' if not torch.cuda.is_available() else 'cuda',
                # freeze bottom layers
                callbacks=[
                    # Freezer(lambda x: not x.startswith('classifier')),
                    EpochScoring("roc_auc", lower_is_better=False,
                                name='val_auroc'),
                    EpochScoring("f1", lower_is_better=False, name='val_f1'),
                    # SacredLogger(,log_on_batch_end=True),
                    LRScheduler(policy='StepLR', step_size=8, gamma=0.2)
                ],
                # shuffle training data
                iterator_train__shuffle=True
            )
            # net.fit(x_train, torch.Tensor(y_train))
            net.set_params(train_split=False, verbose=3)
            if hyperparameter_search_strategy == 'grid':
                grid_search = GridSearchCV(
                    estimator=net, param_grid=params, cv=inner_cv, scoring=('roc_auc'), verbose=3)
            elif hyperparameter_search_strategy == 'random':
                grid_search = RandomizedSearchCV(
                    estimator=net,
                    param_distributions=params,
                    n_iter=random_search_iter,
                    cv=inner_cv,
                    scoring=('roc_auc'),
                    verbose=3,
                    random_state=repeat_seed)
            # grid search using train data
            grid_search.fit(torch.Tensor(x_train.swapaxes(1, 3).swapaxes(2, 3)),
                            torch.Tensor(y_train), groups=z_train)
            print('best params: ', grid_search.best_params_)
            print('best score: ', grid_search.best_score_)
            # # save best params to a file -----------------------------------------------
            best_params_filename = 'repeat_{}_fold_{}_best_params_from_CV.csv'.format(
                which_repeat, which_fold)
            best_params_filename_fullpath = roc_curve_csv_path_ind_test + best_params_filename
            # write best_params to best_params_filename_fullpath
            best_param_wrapped = {k: [v]
                                    for k, v in grid_search.best_params_.items()}
            pd.DataFrame(best_param_wrapped).to_csv(
                best_params_filename_fullpath, index=False)
            print('best params saved to: ', best_params_filename_fullpath)
            # save best score to a file ------------------------------------------------
            best_score_filename = 'repeat_{}_fold_{}_best_score_from_CV.csv'.format(
                which_repeat, which_fold)
            best_score_filename_fullpath = roc_curve_csv_path_ind_test + best_score_filename
            best_score_wrapped = {'best_score': [grid_search.best_score_]}
            pd.DataFrame(best_score_wrapped).to_csv(
                best_score_filename_fullpath, index=False)
            print('best score saved to: ', best_score_filename_fullpath)
            # get the best model -------------------------------------------------------
            net = grid_search.best_estimator_
            # # save the best model to file
            best_model_filename = 'repeat_{}_fold_{}_best_model_from_CV.pkl'.format(
                which_repeat, which_fold)
            best_model_filename_fullpath = roc_curve_csv_path_ind_test + best_model_filename
            # write net to betst_model_filename_fullpath
            with open(best_model_filename_fullpath, 'wb') as f:
                pickle.dump(net, f)
            print('best model saved to: ', best_model_filename_fullpath)

        if not train_from_scratch:
            # count_parameters(model)
            net = NeuralNetBinaryClassifier(
                model,
                max_epochs=model_epochs,
                lr=model_learning_rate,
                batch_size=model_batch_size,
                criterion=nn.BCEWithLogitsLoss(pos_weight=pos_weight),
                optimizer=optim.Adam,
                device='cpu' if not torch.cuda.is_available() else 'cuda',
                # freeze bottom layers
                callbacks=[
                    Freezer(lambda x: not x.startswith('classifier')),
                    EpochScoring("roc_auc"),
                    # SacredLogger(,log_on_batch_end=True),
                    LRScheduler(policy='StepLR', step_size=8, gamma=0.2)
                ],
                # shuffle training data
                iterator_train__shuffle=True
                # class weights
                # criterion__weight=torch.Tensor(list(model_class_weight.values())),
                # validation data
            )
            net.fit(torch.Tensor(x_train.swapaxes(1, 3).swapaxes(2, 3)),
                    torch.Tensor(y_train))
            # print the number of parameters in model
            count_total_parameters(model)
            count_trainable_parameters(model)

            # -------------------------------------------------------------------------
            # fine-tune the model
            # -------------------------------------------------------------------------
            # unfreeze and fine-tune the model
            for p in model.parameters():
                p.requires_grad = True
            count_trainable_parameters(model)
            net.set_params(lr=model_finetune_learning_rate)
            net.set_params(max_epochs=model_fine_tune_epochs)
            net.set_params(callbacks=[])
            net.fit(torch.Tensor(x_train.swapaxes(1, 3).swapaxes(2, 3)),
                    torch.Tensor(y_train))

        ###############################################################################
        # test the model
        ###############################################################################

        y_pred = net.predict_proba(torch.Tensor(
            x_test.swapaxes(1, 3).swapaxes(2, 3)))[:, 1]

        y_pred_class = net.predict(torch.Tensor(
            x_test.swapaxes(1, 3).swapaxes(2, 3)))

        ###############################################################################
        # metrics
        ###############################################################################

        tn, fp, fn, tp = confusion_matrix(y_test, y_pred_class).ravel()
        print(confusion_matrix(y_test, y_pred_class))
        print('tn: ', tn)
        print('fp: ', fp)
        print('fn: ', fn)
        print('tp: ', tp)

        # print classification report
        print(classification_report(y_test, y_pred_class))

        # make a dataframe including  cols, y_test, y_pred, z_test
        y_pred_f = y_pred.flatten()
        y_pred_class_f = y_pred_class.flatten()
        y_test_f = y_test.flatten()
        z_test_f = z_test.flatten()
        bam_id_test_f = bam_id_test.flatten()

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

        # if_set_sample_weight_for_testing = True, calculate the sample_weight for y_test
        if if_set_sample_weight_for_testing:
            from sklearn.utils import class_weight
            sample_weight = class_weight.compute_sample_weight(
                class_weight='balanced', y=y_test)
        else:
            sample_weight = None

        # add repeat number
        results['repeat'] = which_repeat
        # add fold number
        results['fold'] = which_fold

        # add model_name
        results['model_name'] = model_name
        # calculate the the auroc
        auroc = roc_auc_score(y_test, y_pred, sample_weight=sample_weight)
        results['test_auroc'] = auroc
        print('auroc: ', auroc)

        # calculate the the auprc
        auprc = average_precision_score(
            y_test, y_pred, sample_weight=sample_weight)
        results['test_auprc'] = auprc
        print('auprc: ', auprc)

        # calculate generic accuracy
        acc = accuracy_score(y_test, y_pred > 0.5, sample_weight=sample_weight)
        results['test_acc'] = acc
        print('acc: ', acc)

        # calculate generic ppv
        ppv = precision_score(
            y_test, y_pred > 0.5, sample_weight=sample_weight, zero_division=0)
        results['test_ppv'] = ppv
        print('ppv: ', ppv)

        # calculate generic sensitivity
        sensitivity = recall_score(
            y_test, y_pred > 0.5, sample_weight=sample_weight)
        results['test_sensitivity'] = sensitivity
        print('sensitivity: ', sensitivity)

        # calculate generic specificity
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
        roc_filename = 'repeat_{}_fold_{}_roc.csv'.format(
            which_repeat, which_fold)
        roc_filename_fullpath = roc_curve_csv_path_ind_test + roc_filename
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

        tn_sample_namelist = bam_id_test[tn_sample_id]
        fp_sample_namelist = bam_id_test[fp_sample_id]
        fn_sample_namelist = bam_id_test[fn_sample_id]
        tp_sample_namelist = bam_id_test[tp_sample_id]

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

        # -------------------------------------------------------------------------------
        # get the sensitivity at 95% specificity, and threshold
        # -------------------------------------------------------------------------------
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

# -------------------------------------------------------------------------------
# Append values to the lists associated with each key
# -------------------------------------------------------------------------------
        for key, value in results.items():
            results_container[key].append(value)
        # convert results to a dictionary
        results_dict = dict(results_container)
        results_df = pd.DataFrame.from_dict(results_dict)
        # convet to dataframe
        metrics_filename = 'repeat_{}_fold_{}_metrics.csv'.format(
            which_repeat, which_fold)
        metrics_filename_fullpath = roc_curve_csv_path_ind_test + metrics_filename
        results_df.to_csv(metrics_filename_fullpath, index=False)
        print('metrics saved to: ', metrics_filename_fullpath)
        # print key-value pairs in results
        for key, value in results.items():
            print(key, value)

        print('* ' * 50)
