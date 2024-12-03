import cupy as cp
from matplotlib import pyplot as plt
import seaborn as sns
import shap
import pickle
import random
import argparse
import os
from collections import defaultdict
import sklearn.model_selection
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.model_selection import KFold, StratifiedKFold, StratifiedGroupKFold

from sklearn.metrics import roc_curve, precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score

from sklearn.utils import class_weight

from sklearn.datasets import make_classification
from numpy import std
from numpy import mean
import numpy as np
import xgboost as xgb
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
# shap

# gpu
################################################################################
# parameters
################################################################################
parser = argparse.ArgumentParser(
    description='XGBoost Binary Function Nested CV')
parser.add_argument('--num_repeats', type=int,
                    default=10, help='Number of repeats')
parser.add_argument('--outer_cv_fold', type=int, default=5,
                    help='Number of outer cross-validation folds')
parser.add_argument('--inner_cv_fold', type=int, default=3,
                    help='Number of inner cross-validation folds')
parser.add_argument('--seed', type=int, default=0, help='Random seed')
parser.add_argument('--xgb_training_device', type=str,
                    default="cuda", help='Training device ("cuda" or "cpu")')
parser.add_argument('--if_set_sample_weight_for_testing', type=bool,
                    default=False, help='Whether to set sample weight for testing')
parser.add_argument('--input_csv', type=str, default="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb9/0-0.03_pairwise/xgboost_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0-0.03_pairwise.csv", help='Input CSV file path')

# add param objective='binary:logistic'
parser.add_argument('--xgb_objective', type=str, default='binary:logistic')
# tree_method='hist'
parser.add_argument('--xgb_tree_method', type=str, default='hist')
# nthread=-1
parser.add_argument('--xgb_nthread', type=int, default=-1)
# add param shap_max_display=30
parser.add_argument('--shap_max_display', type=int, default=30,
                    help='Maximum number of features to display in SHAP summary plot')

# add a param called "quick_test_run", with default to False
parser.add_argument('--quick_test_run', action='store_true',
                    help='Whether to run quick test')

# add a param called "hyperparameter_search_strategy", with default to "grid", choices=["grid", "random"]
parser.add_argument('--hyperparameter_search_strategy', type=str, default='grid',
                    choices=['grid', 'random'], help='Hyperparameter search strategy')
# add a param called "random_search_iter", with default to 100
parser.add_argument('--random_search_iter', type=int,
                    default=100, help='Number of iterations for random search')

parser.add_argument('--keep_cnv_19', action='store_true')


args = parser.parse_args()

seed = args.seed
if_set_sample_weight_for_testing = args.if_set_sample_weight_for_testing
input_csv = args.input_csv
keep_cnv_19 = args.keep_cnv_19


# cross-validation parameters
num_repeats = args.num_repeats
outer_cv_fold = args.outer_cv_fold
inner_cv_fold = args.inner_cv_fold
quick_test_run = args.quick_test_run
# hyperparameter search strategy
hyperparameter_search_strategy = args.hyperparameter_search_strategy
random_search_iter = args.random_search_iter
# model parameters
xgb_training_device = args.xgb_training_device
xgb_objective = args.xgb_objective
xgb_tree_method = args.xgb_tree_method
xgb_nthread = args.xgb_nthread

# shap parameters
shap_max_display = args.shap_max_display


################################################################################
# create a folder to store the roc curve csv
basename = os.path.basename(input_csv)
basename2 = os.path.splitext(basename)[0]
roc_curve_csv_path = os.path.dirname(
    input_csv) + '/' + basename2 + '_roc_curve_csv/'
os.makedirs(roc_curve_csv_path, exist_ok=True)

################################################################################
print('read input csv...')
# read input csv
data = pd.read_csv(input_csv)
# summarize the values in 'bicohort' column


# label the bicohort column, Healthy is 0, Cancer is 1
data['bicohort'] = data['bicohort'].replace(
    'Healthy', 0).infer_objects(copy=False)
data['bicohort'] = data['bicohort'].replace(
    'Cancer', 1).infer_objects(copy=False)

colnames = data.columns

# get the x, y
# y is the bicohort column
# x is the rest of the columns
# if the column is 'primary', 'bicohort', 'patient_id', drop it
x = data.drop(['primary', 'bicohort', 'patient_id'], axis=1)
y = data['bicohort']
z = data['patient_id']
bam_id = data['primary']


# filter out the columns names containing 'cnv.19' unless keep_cnv_19 is True
if not keep_cnv_19:
    x = x.loc[:, ~x.columns.str.contains('cnv.19')]
    x = x.loc[:, ~x.columns.str.contains('ctRatio.19')]
    x = x.loc[:, ~x.columns.str.contains('slRatio.19')]
for feature in x.columns:
    print('feature:', feature)


model = xgb.XGBClassifier(
    objective=xgb_objective,
    tree_method=xgb_tree_method,
    nthread=xgb_nthread,
    device=xgb_training_device,
    seed=seed)

# https://randomrealizations.com/posts/xgboost-parameter-tuning-with-optuna/
xgb_params = {
    # boosting parameters
    'learning_rate': [0.01, 0.001],
    'n_estimators': [500, 1000],
    # 'early_stopping_rounds': [50],

    # sampling parameters
    'subsample': [0.7, 0.8, 1.0],
    'colsample_bylevel': [0.6, 0.8, 1.0],

    # regularization parameters to prevent overfitting
    # 'reg_lambda': [0, 0.1, 1],

    # tree complexity parameters
    'max_depth': [3, 5],
    'min_child_weight': [1, 5, 10]
}

if quick_test_run:
    xgb_params = {
        # boosting parameters
        'learning_rate': [0.01, 0.001],
        'n_estimators': [500, 1000],
        'max_depth': [3, 5],
        'min_child_weight': [1, 5, 10]
    }


# Initialize lists to store the evaluation results
outer_scores = []
inner_scores = []
results = {}
results_container = defaultdict(list)

SHAP_values_per_fold = []

shap_values_per_cv = dict()
for sample in x.index:
    # Create keys for each sample
    shap_values_per_cv[sample] = {}
    # Then, keys for each CV fold within each sample
    for CV_repeat in range(num_repeats):
        shap_values_per_cv[sample][CV_repeat+1] = {}


shap_values_per_cv_with_bam = dict()
for bam in bam_id:
    # Create keys for each sample
    shap_values_per_cv_with_bam[bam] = {}
    # Then, keys for each CV fold within each sample
    for CV_repeat in range(num_repeats):
        shap_values_per_cv_with_bam[bam][CV_repeat+1] = {}


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


################################################################################
# nested cv start
################################################################################
np.random.seed(seed)

# print the parameters
print('num_repeats: ', num_repeats)
print('outer_cv_fold: ', outer_cv_fold)
print('inner_cv_fold: ', inner_cv_fold)
print('seed: ', seed)
print('xgb_training_device: ', xgb_training_device)
print('xgb_objective: ', xgb_objective)
print('xgb_tree_method: ', xgb_tree_method)
print('xgb_nthread: ', xgb_nthread)
print('shap_max_display: ', shap_max_display)
print('if_set_sample_weight_for_testing: ', if_set_sample_weight_for_testing)
print('input_csv: ', input_csv)
print('roc_curve_csv_path: ', roc_curve_csv_path)
print('x shape: ', x.shape)
print('y shape: ', y.shape)
print('z shape: ', z.shape)
print('bam_id shape: ', bam_id.shape)
print('start the repeats...')

for i in range(num_repeats):
    # Generate a new seed for each repeat
    repeat_seed = seed + i
    np.random.seed(repeat_seed)
    # tf.random.set_seed(repeat_seed)
    random.seed(repeat_seed)

    # report which loop we are on
    which_repeat = i + 1
    print('= ' * 50)
    print('Repeat: ', which_repeat)
    print('Seed: ', repeat_seed)

    outer_cv = StratifiedGroupKFold(
        n_splits=outer_cv_fold, shuffle=True, random_state=repeat_seed)

    # split into train and test
    for fold_i, (train_idx, test_idx) in enumerate(outer_cv.split(x, y, groups=z)):
        # print repeat and fold number
        which_fold = fold_i + 1
        print('- ' * 50)
        print('Repeat: ', which_repeat)
        print('Fold: ', which_fold)
        print('Seed: ', repeat_seed)
        # Split the data into train and test sets
        x_train, x_test = x.iloc[train_idx], x.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        z_train, z_test = z.iloc[train_idx], z.iloc[test_idx]
        bam_id_train, bam_id_test = bam_id.iloc[train_idx], bam_id.iloc[test_idx]
        within_fold_classes_weights = class_weight.compute_sample_weight(
            class_weight='balanced',
            y=y_train
        )

        # convert into cupy array
        if xgb_training_device == 'cuda':
            x_train_cp = cp.array(x_train.to_numpy()).get()
            x_test_cp = cp.array(x_test.to_numpy()).get()
            y_train_cp = cp.array(y_train.to_numpy()).get()
            y_test_cp = cp.array(y_test.to_numpy()).get()
        elif xgb_training_device == 'cpu':
            x_train_cp = x_train
            x_test_cp = x_test
            y_train_cp = y_train
            y_test_cp = y_test
        else:
            raise ValueError(
                'xgb_training_device should be either "cuda" or "cpu"')

        # prepare inner cross-validation, repeat 5 times
        inner_cv = StratifiedGroupKFold(
            n_splits=inner_cv_fold, shuffle=True, random_state=repeat_seed)

        # Perform grid search for hyperparameter tuning if hyperparameter_search_strategy is "grid"
        if hyperparameter_search_strategy == 'grid':
            grid_search = GridSearchCV(
                estimator=model, param_grid=xgb_params, cv=inner_cv, scoring=('roc_auc'), verbose=3)
        elif hyperparameter_search_strategy == 'random':
            grid_search = RandomizedSearchCV(
                estimator=model, param_distributions=xgb_params, n_iter=random_search_iter, cv=inner_cv, scoring=('roc_auc'), verbose=3, random_state=repeat_seed)

        grid_search.fit(x_train_cp, y_train_cp, groups=z_train,
                        sample_weight=within_fold_classes_weights)

        # save the grid_search object to a pickle file
        grid_search_pickle_path = roc_curve_csv_path + \
            'repeat_{}_fold_{}_grid_search.pkl'.format(
                which_repeat, which_fold)
        with open(grid_search_pickle_path, 'wb') as f:
            pickle.dump(grid_search, f)

        ########################################################################
        # get the best mode and predict
        ########################################################################
        # Get the best model from grid search
        best_model = grid_search.best_estimator_
        # using test data
        print("Evaluate...")
        # model prediction
        y_pred_class = best_model.predict(x_test_cp)
        y_pred = best_model.predict_proba(x_test_cp)[:, 1]
        ########################################################################
        # SHAP
        ########################################################################
        # calculate and store SHAP value
        explainer = shap.TreeExplainer(best_model)
        shap_values = explainer.shap_values(x_test_cp)

        shap_values_filename = 'repeat_{}_fold_{}_SHAP.csv'.format(
            which_repeat, which_fold)
        shap_values_filename_fullpath = roc_curve_csv_path + shap_values_filename
        shap_values_df = pd.DataFrame(shap_values, columns=x_test.columns)
        shap_values_df.to_csv(shap_values_filename_fullpath, index=False)
        print('fold shap values saved to: ', shap_values_filename_fullpath)

        # -----------------------------------------------------------------------
        # plot the shap summary plot, dot plot
        # -----------------------------------------------------------------------
        shap.summary_plot(np.array(shap_values), x_test,
                          max_display=shap_max_display, show=False, plot_type='dot')
        # make a path to store the shap summary plot
        within_fold_shap_plot_path = roc_curve_csv_path + \
            'repeat_{}_fold_{}_SHAP_summary_plot_dot.pdf'.format(
                which_repeat, which_fold)
        plot_title = 'Repeat: {}, Fold: {}'.format(which_repeat, which_fold)
        plt.title(plot_title)
        plt.savefig(within_fold_shap_plot_path)
        plt.clf()
        message = 'fold shap summary plot saved to: {}'.format(
            within_fold_shap_plot_path)
        print(message)

        # -----------------------------------------------------------------------
        # plot the shap summary plot, bar plot
        # -----------------------------------------------------------------------
        shap.summary_plot(np.array(shap_values), x_test,
                          max_display=shap_max_display, show=False, plot_type='bar')
        # make a path to store the shap summary plot
        within_fold_shap_plot_path = roc_curve_csv_path + \
            'repeat_{}_fold_{}_SHAP_summary_plot_bar.pdf'.format(
                which_repeat, which_fold)
        plot_title = 'Repeat: {}, Fold: {}'.format(which_repeat, which_fold)
        plt.title(plot_title)
        plt.savefig(within_fold_shap_plot_path)
        plt.clf()
        message = 'fold shap summary plot saved to: {}'.format(
            within_fold_shap_plot_path)
        print(message)

        # -----------------------------------------------------------------------
        # plot the shap summary plot, violin plot - removed due to dimension issue in small sample groups
        # -----------------------------------------------------------------------
        # shap.summary_plot(np.array(shap_values), x_test,
        #                   max_display=shap_max_display, show=False, plot_type='violin')
        # # make a path to store the shap summary plot
        # within_fold_shap_plot_path = roc_curve_csv_path + \
        #     'repeat_{}_fold_{}_SHAP_summary_plot_violin.pdf'.format(
        #         which_repeat, which_fold)
        # plot_title = 'Repeat: {}, Fold: {}'.format(which_repeat, which_fold)
        # plt.title(plot_title)
        # plt.savefig(within_fold_shap_plot_path)
        # plt.clf()
        # message = 'fold shap summary plot saved to: {}'.format(
        #     within_fold_shap_plot_path)
        # print(message)

        # store the shap values in all repeats

        for i, sample_id in enumerate(test_idx):
            shap_values_per_cv[sample_id][which_repeat] = shap_values[i]

        for i, which_bam_id in enumerate(bam_id_test):
            # which_bam = bam_id_test.iloc[which_bam_id]
            shap_values_per_cv_with_bam[which_bam_id][which_repeat] = shap_values[i]

        ########################################################################
        # metrics
        ########################################################################
        from sklearn.metrics import confusion_matrix
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred_class).ravel()
        print(confusion_matrix(y_test, y_pred_class))
        print('tn: ', tn)
        print('fp: ', fp)
        print('fn: ', fn)
        print('tp: ', tp)

        # print classification report
        from sklearn.metrics import classification_report
        print(classification_report(y_test, y_pred_class))

        # make a dataframe including  cols, y_test, y_pred, z_test
        y_pred_f = y_pred.flatten()
        y_pred_class_f = y_pred_class.flatten()
        y_test_f = y_test.to_numpy().flatten()
        z_test_f = z_test.to_numpy().flatten()
        bam_id_test_f = bam_id_test.to_numpy().flatten()

        # create a dataframe to store y_test_f, y_pred_f, z_test_f
        y_pred_df = pd.DataFrame(
            {'repeat': which_repeat,
             'fold': which_fold,
             'model_name': 'xgboost',
             'y_test': y_test_f,
             'y_pred': y_pred_f,
             'y_pred_class': y_pred_class_f,
             'z_test': z_test_f,
             'bam_id_test': bam_id_test_f})

        y_pred_df_filename = 'repeat_{}_fold_{}_y_test_y_pred_z_test.csv'.format(
            which_repeat, which_fold)
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
        results['repeat'] = which_repeat
        # add fold number
        results['fold'] = which_fold

        # add model name
        results['model_name'] = 'xgboost'

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
            y_test, y_pred > 0.5, 
            sample_weight=sample_weight,
            zero_division=0)
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
        roc_filename = 'repeat_{}_fold_{}_roc.csv'.format(
            which_repeat, which_fold)
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
            which_repeat, which_fold)
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

        tn_sample_namelist = bam_id_test.iloc[tn_sample_id]
        fp_sample_namelist = bam_id_test.iloc[fp_sample_id]
        fn_sample_namelist = bam_id_test.iloc[fn_sample_id]
        tp_sample_namelist = bam_id_test.iloc[tp_sample_id]

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
            which_repeat, which_fold)
        metrics_filename_fullpath = roc_curve_csv_path + metrics_filename
        results_df.to_csv(metrics_filename_fullpath, index=False)
        print('metrics saved to: ', metrics_filename_fullpath)

        print(results_df)

################################################################################
# handle shap values
################################################################################

average_shap_values, bam_id_col, stds, ranges = [], [], [], []

for bam in bam_id:
    # Get all SHAP values for sample number i
    df_per_obs = pd.DataFrame.from_dict(shap_values_per_cv_with_bam[bam])
    # Get relevant statistics for every sample
    average_shap_values.append(df_per_obs.mean(axis=1).values)
    # add bam_id to average_shap_values as the first column
    bam_id_col.append([bam] * len(df_per_obs.mean(axis=1).values))
    stds.append(df_per_obs.std(axis=1).values)
    ranges.append(df_per_obs.max(axis=1).values-df_per_obs.min(axis=1).values)


# convert average_shap_values, bam_id_col, stds, ranges to dictionary, save to pickle file

shap_values_dict = {'shap_values_per_cv_with_bam': shap_values_per_cv_with_bam,
                    'average_shap_values': average_shap_values,
                    'bam_id_col': bam_id_col,
                    'stds': stds,
                    'ranges': ranges}

shap_values_dict_filename = 'SHAP.all_repeat_shap_values_dict.pkl'
shap_values_dict_filename_fullpath = roc_curve_csv_path + shap_values_dict_filename
with open(shap_values_dict_filename_fullpath, 'wb') as f:
    pickle.dump(shap_values_dict, f)
print('shap values saved to: ', shap_values_dict_filename_fullpath)

# plot

shap.summary_plot(np.array(average_shap_values), x,
                  max_display=30, show=False, plot_type='dot')
title = 'Average SHAP values after cross-validation, N repeat:{}, N fold:{}'.format(
    num_repeats, outer_cv_fold)
plt.title(title)
average_SHAP_plot_path = roc_curve_csv_path + \
    'SHAP.average_SHAP_summary_plot.dot.pdf'
plt.savefig(average_SHAP_plot_path)
plt.clf()
print('Average SHAP dot plot saved to: ', average_SHAP_plot_path)

# -------------------------------------------------------------------------------
# bar plot
# -------------------------------------------------------------------------------
shap.summary_plot(np.array(average_shap_values), x,
                  max_display=30, show=False, plot_type='bar')
title = 'Average SHAP values after cross-validation, N repeat:{}, N fold:{}'.format(
    num_repeats, outer_cv_fold)
plt.title(title)
average_SHAP_plot_path = roc_curve_csv_path + \
    'SHAP.average_SHAP_summary_plot.bar.pdf'
plt.savefig(average_SHAP_plot_path)
plt.clf()
print('Average SHAP bar plot saved to: ', average_SHAP_plot_path)


################################################################################
print('Jobs done!')
################################################################################
