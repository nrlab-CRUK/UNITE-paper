from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.datasets import make_classification
from numpy import std
from numpy import mean
import numpy as np
import xgboost as xgb
from sklearn.model_selection import GridSearchCV, StratifiedKFold
import pandas as pd
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
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
import os
import argparse
import random
import pickle
# shap
import shap
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
import cupy as cp
################################################################################
# parameters
################################################################################
parser = argparse.ArgumentParser(
    description='XGBoost Binary Function Nested CV')
parser.add_argument('--num_repeats', type=int,
                    default=20, help='Number of repeats')
parser.add_argument('--outer_cv_fold', type=int, default=5,
                    help='Number of outer cross-validation folds')
parser.add_argument('--inner_cv_fold', type=int, default=3,
                    help='Number of inner cross-validation folds')
parser.add_argument('--seed', type=int, default=0, help='Random seed')
parser.add_argument('--xgb_training_device', type=str,
                    default="cuda", help='Training device ("cuda" or "cpu")')
parser.add_argument('--if_set_sample_weight_for_testing', type=bool,
                    default=False, help='Whether to set sample weight for testing')
parser.add_argument('--input_csv', type=str, default="/scratchc/nrlab/wang04/ulyses_results_iteration/iteration_3_merge_below_3/xgboost_model/xgb6/0-0.03_pairwise/xgboost_input.feat.cnv_sd_length_ctRatio_slRatio.tf.0-0.03_pairwise.csv", help='Input CSV file path')

# add param objective='binary:logistic'
parser.add_argument('--xgb_objective', type=str, default='binary:logistic')
# tree_method='hist'
parser.add_argument('--xgb_tree_method', type=str, default='hist')
# nthread=-1
parser.add_argument('--xgb_nthread', type=int, default=-1)


args = parser.parse_args()

seed = args.seed
if_set_sample_weight_for_testing = args.if_set_sample_weight_for_testing
input_csv = args.input_csv
# cross-validation parameters
num_repeats = args.num_repeats
outer_cv_fold = args.outer_cv_fold
inner_cv_fold = args.inner_cv_fold
# model parameters
xgb_training_device = args.xgb_training_device
xgb_objective = args.xgb_objective
xgb_tree_method = args.xgb_tree_method
xgb_nthread = args.xgb_nthread


################################################################################
# create a folder to store the roc curve csv
basename = os.path.basename(input_csv)
basename2 = os.path.splitext(basename)[0]
roc_curve_csv_path = os.path.dirname(
    input_csv) + '/' + basename2 + '_roc_curve_csv/'
os.makedirs(roc_curve_csv_path, exist_ok=True)

################################################################################
# read input csv
data = pd.read_csv(input_csv)
# summarize the values in 'bicohort' column


# label the bicohort column, Healthy is 0, Cancer is 1
data['bicohort'] = data['bicohort'].replace('Healthy', 0)
data['bicohort'] = data['bicohort'].replace('Cancer', 1)

colnames = data.columns

# get the x, y
# y is the bicohort column
# x is the rest of the columns
# if the column is 'primary', 'bicohort', 'patient_id', drop it
x = data.drop(['primary', 'bicohort', 'patient_id'], axis=1)
y = data['bicohort']
z = data['patient_id']
bam_id = data['primary']

classes_weights = dict(
    zip(np.unique(y),
        class_weight.compute_class_weight(class_weight='balanced',
                                          classes=np.unique(np.array(y)), y=y)))
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


    # sampling parameters
    'subsample': [0.7, 0.8, 1.0],
    'colsample_bylevel': [0.6, 0.8, 1.0],

    # regularization parameters to prevent overfitting
    # 'reg_lambda': [0, 0.1, 1],

    # tree complexity parameters
    'max_depth': [3, 5],
    'min_child_weight': [1, 5, 10]
}


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


xgb_training_device = 'cuda'
model = xgb.XGBClassifier(
    objective=xgb_objective,
    tree_method=xgb_tree_method,
    nthread=xgb_nthread,
    device=xgb_training_device,
    seed=seed)


repeat_seed = 1


# split x into 70% train and 30% validation, get the split index


x_train, x_test, y_train, y_test, z_train, z_test, bam_id_train, bam_id_test = train_test_split(
    x, y, z, bam_id, test_size=0.3, random_state=42)

# convert x_train, x_test, y_train, y_test to numpy array


x_train_cp = cp.array(x_train.to_numpy()).get()
x_test_cp = cp.array(x_test.to_numpy()).get()
y_train_cp = cp.array(y_train.to_numpy()).get()
y_test_cp = cp.array(y_test.to_numpy()).get()

# prepare inner cross-validation, repeat 5 times
inner_cv = StratifiedGroupKFold(
    n_splits=inner_cv_fold, shuffle=True, random_state=repeat_seed)

# Perform grid search for hyperparameter tuning
grid_search = GridSearchCV(
    estimator=model, param_grid=xgb_params, cv=inner_cv, scoring=('roc_auc'), verbose=3)

classes_weights = class_weight.compute_sample_weight(
    class_weight='balanced',
    y=y_train
)

grid_search.fit(x_train_cp, y_train_cp, groups=z_train,
                sample_weight=classes_weights)

# save the grid_search object to a pickle file
i = 0
fold_i = 1
grid_search_pickle_path = roc_curve_csv_path + \
    '/scratchc/nrlab/wang04/ulyses/repeat_{}_fold_{}_grid_search.pkl'.format(
        i + 1, fold_i + 1)

grid_search_pickle_path = '/home/nrlab/wang04/ulyses/repeat_{}_fold_{}_grid_search.pkl'.format(
    i + 1, fold_i + 1)

with open(grid_search_pickle_path, 'wb') as f:
    pickle.dump(grid_search, f)

best_model = grid_search.best_estimator_
# using test data
print("Evaluate...")
# model prediction
y_pred = best_model.predict(x_test_cp)
########################################################################
# SHAP
########################################################################
# calculate and store SHAP value
# explainer = shap.Explainer(best_model)
explainer = shap.TreeExplainer(best_model)
shap_values = explainer.shap_values(x_test_cp)
# for SHAPs in shap_values:
#     SHAP_values_per_fold.append(SHAPs)

shap_values_filename = 'repeat_{}_fold_{}_SHAP.csv'.format(
    i + 1, fold_i + 1)
shap_values_filename_fullpath = roc_curve_csv_path + shap_values_filename
shap_values_df = pd.DataFrame(shap_values, columns=x_test.columns)
shap_values_df.to_csv(shap_values_filename_fullpath, index=False)
print('fold shap values saved to: ', shap_values_filename_fullpath)


# plt.rcParams['font.size'] = 5
shap_max_display = 30
shap.summary_plot(np.array(shap_values), x_test,
                  max_display=shap_max_display, show=False, plot_type='dot')
# plt.title('Average SHAP values after cross-validation')
# set the plot size: width = 90mm, height = 200mm
# plt.gcf().set_size_inches(w=90/25.4, h=200/25.4)
plt.savefig('/home/nrlab/wang04/ulyses/shap_summary_plot.pdf')
# clear the plot
plt.clf()

shap.summary_plot(np.array(shap_values), x_test,
                  max_display=shap_max_display, show=False, plot_type='bar')
plt.savefig('/home/nrlab/wang04/ulyses/shap_summary_plot_bar.pdf')
plt.clf()


shap.summary_plot(np.array(shap_values), x_test,
                  max_display=shap_max_display, show=False, plot_type='dot')
plt.savefig('/home/nrlab/wang04/ulyses/shap_summary_plot_dot.pdf')
plt.clf()


shap.summary_plot(np.array(shap_values), x_test,
                  max_display=shap_max_display, show=False, plot_type='violin')
plt.savefig('/home/nrlab/wang04/ulyses/shap_summary_plot_violin.pdf')
plt.clf()

# independent instance plot

# shap waterfall plot
shap_values_obj = explainer(x_test)
shap.plots.waterfall(shap_values_obj[0], show=False)
plt.savefig('/home/nrlab/wang04/ulyses/shap_waterfall_plot_onesample.pdf')
plt.clf()

shap.plots.bar(shap_values_obj[0], show=False)
plt.savefig('/home/nrlab/wang04/ulyses/shap_bar_plot_onesample.pdf')
plt.clf()
