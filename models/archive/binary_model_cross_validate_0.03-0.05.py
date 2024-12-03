

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import sklearn as sk


import numpy as np
import matplotlib as mpl


# cv parameters
num_folds = 5 
num_repeats = 10 

# model parameters
optimizer = 'adam'
learning_rate = 0.0008
loss = 'binary_crossentropy'

# training parameters
epochs = 35 
batch_size = 16 
metrics=['accuracy',
        tf.keras.metrics.AUC(curve = "ROC", name='AUROC'),
        tf.keras.metrics.AUC(curve = "PR", name='AUPRC'),
         tf.keras.metrics.Precision(name='precision'),
         tf.keras.metrics.Recall(name='recall'),
         tf.keras.metrics.SensitivityAtSpecificity(specificity=0.95, name="sen_at_0.95spec"),
         tf.keras.metrics.SensitivityAtSpecificity(specificity=0.98, name="sen_at_0.98spec")]

# readin data

x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.03_0.05_.npy'

x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.01_0.03_.npy'

x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0_0.01_.npy'


x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.05_0.1_.npy'

x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.1_0.2_.npy'

x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.npy'


y_file = 'y' + x_file[1:]
z_file = 'z' + x_file[1:]
performance_csv = 'performance' + x_file[7:-4] + '.debug.csv'

# load the data
x = np.load(x_file, allow_pickle=True)
y = np.load(y_file, allow_pickle=True)
z = np.load(z_file, allow_pickle=True)

y = np.array(y)
input_shape = x.shape[1:]

# define model
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier
import scikeras.wrappers as sw


def create_model(input_shape, loss, learning_rate, metrics):
    model = keras.Sequential(
        [
            keras.Input(shape=input_shape),
            
            layers.Conv2D(16, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(3, 2)),
            
            layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),
            
            layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),
            
            layers.Conv2D(64, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),
            
            layers.Flatten(),
            layers.Dense(32, activation="relu"),
            #layers.Dropout(0.2),
            layers.Dense(1, activation="sigmoid")
        ])
    
    opt = keras.optimizers.Adam(learning_rate=learning_rate)
    model.compile(loss=loss, 
        optimizer=opt, 
        metrics=metrics)
    return model


# compute class weights and set dictionary to class_weight parameter
from sklearn.utils import class_weight
class_weight = dict(zip(np.unique(y), class_weight.compute_class_weight(class_weight = 'balanced', 
                                                                        classes = np.unique(y), 
                                                                        y = y)))

model = sw.KerasClassifier(
    model=create_model,
    model__metrics=metrics,
    model__input_shape=input_shape,
    model__loss=loss,
    model__learning_rate=learning_rate,
    class_weight=class_weight,
    validation_split=0.2,
    epochs=epochs,
    batch_size=batch_size,
    verbose=1
)


# define scores
from sklearn.metrics import make_scorer
from sklearn.metrics import recall_score
from sklearn.metrics import confusion_matrix

#TODO this is not working
def sensitivity_at_specificity(y_true, y_pred_proba, specificity):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_proba >= 0.5).ravel()
    actual_negatives = tn + fp
    actual_positives = fn + tp
    tn_rate = tn / actual_negatives
    fp_rate = fp / actual_negatives
    cutoff = 0.5
    while fp_rate > (1-specificity) and cutoff > 0:
        cutoff -= 0.01
        y_pred = (y_pred_proba > cutoff).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        actual_negatives = tn + fp
        actual_positives = fn + tp
        tn_rate = tn / actual_negatives
        fp_rate = fp / actual_negatives
    sensitivity = tp / actual_positives
    return sensitivity

# define a function to calculate the sensitivity at a given specificity
def sensitivity_at_specificity(y_true, y_pred_proba, specificity):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred_proba >= 0.5).ravel()
    actual_negatives = tn + fp
    actual_positives = fn + tp
    tn_rate = tn / actual_negatives
    fp_rate = fp / actual_negatives
    cutoff = 0.5
    while fp_rate > (1-specificity) and cutoff > 0:
        cutoff -= 0.01
        y_pred = (y_pred_proba > cutoff).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        actual_negatives = tn + fp
        actual_positives = fn + tp
        tn_rate = tn / actual_negatives
        fp_rate = fp / actual_negatives
    sensitivity = tp / actual_positives
    return sensitivity


scoring = {'specificity': make_scorer(recall_score, pos_label=0),
	'sensitivity': make_scorer(recall_score, pos_label=1),
	'acc': 'accuracy',
	'f1': 'f1',
	'ppv': 'precision',
	'roc_auc': 'roc_auc',
	'pr_auc': 'average_precision',
	'sen_at_0.95spec': make_scorer(sensitivity_at_specificity,  specificity=0.95, needs_proba=True),
	'sen_at_0.98spec': make_scorer(sensitivity_at_specificity, specificity=0.98, needs_proba=True)}

# cross validate


from sklearn.model_selection import RepeatedStratifiedKFold 
from sklearn.model_selection import cross_validate
from sklearn.model_selection import StratifiedGroupKFold

kfold_grouped = StratifiedGroupKFold(n_splits=num_folds,
                            shuffle=True)

results = {}

for i in range(num_repeats):
    # report which loop we are on
    print('Repeat: ', i+1)
    
    # Perform cross-validation with shuffling enabled
    cv_results = cross_validate(model, 
                                x, 
                                y, 
                                groups=z, 
                                return_train_score=True,
                                scoring=scoring, 
                                cv=kfold_grouped, 
                                n_jobs=1, 
                                verbose=0)
    
    # Store the performance metrics for each repetition, store both train and test results
    for key, value in cv_results.items():
        if key not in results:
            results[key] = []
        results[key].append(value)

# Create a pandas DataFrame from the results dictionary
import pandas as pd

# remove 'estimator' key
relevant_keys = [k for k in results.keys() if k != 'estimator']
filtered_dict = {k: np.concatenate(v) for k, v in results.items() if k in relevant_keys}

# create a pandas DataFrame from the new dictionary
cv_df = pd.DataFrame.from_dict(filtered_dict)
# save to csv
cv_df.to_csv(performance_csv, index=False)

# clean up memory
del cv_results
del results
tf.keras.backend.clear_session()


