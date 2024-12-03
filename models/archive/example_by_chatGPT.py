



import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import sklearn as sk
import numpy as np

x_file = 'x_clean_n_isize_n_motif_s1_C_n_motif_s1_T_TF_0.2_1_.npy'
y_file = 'y' + x_file[1:]
z_file = 'z' + x_file[1:]
performance_csv = 'performance' + x_file[7:-4] + '.debug.csv'

# load the data
x = np.load(x_file, allow_pickle=True)
y = np.load(y_file, allow_pickle=True)
z = np.load(z_file, allow_pickle=True)
y = np.array(y)

# cv parameters
num_folds = 5 
num_repeats = 2 

# model parameters
model_input_shape = x.shape[1:]
model_learning_rate = 0.0008
model_loss = 'binary_crossentropy'

# training parameters
model_epochs = 35
model_batch_size = 16


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


# define function to calculate metrics
import keras.backend as K

def sens(y_true, y_pred): 
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())

def spec(y_true, y_pred):
    true_negatives = K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1)))
    possible_negatives = K.sum(K.round(K.clip(1 - y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())


model_metrics=[
        'accuracy',
        tf.keras.metrics.AUC(curve = "ROC", name='AUROC'),
        tf.keras.metrics.AUC(curve = "PR", name='AUPRC'),
        tf.keras.metrics.Precision(name='PPV'),
        tf.keras.metrics.Recall(name='recall_sens'),
        sens,
        spec,
        tf.keras.metrics.SensitivityAtSpecificity(specificity=0.95, name="sen_95spe"),
        tf.keras.metrics.SensitivityAtSpecificity(specificity=0.98, name="sen_98spe"),
        tf.keras.metrics.SensitivityAtSpecificity(specificity=0.99, name="sen_99spe")
        ]
# compute class weights and set dictionary to class_weight parameter
from sklearn.utils import class_weight
class_weight = dict(zip(np.unique(y), class_weight.compute_class_weight(class_weight = 'balanced', 
                                                                        classes = np.unique(y), 
                                                                        y = y)))

from sklearn.model_selection import StratifiedGroupKFold
kfold_grouped = StratifiedGroupKFold(n_splits=num_folds,
                            shuffle=True)

results = {}

for i in range(num_repeats):
    # report which loop we are on
    print('#'*50)
    print('Repeat: ', i+1)
    
    # split into train and test
    for train_idx, test_idx in kfold_grouped.split(x, y, groups= z):
        # Split the data into train and test sets
        X_train, X_test = x[train_idx], x[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Create and train the model
        model = create_model(input_shape=model_input_shape, 
                            loss=model_loss, 
                            learning_rate=model_learning_rate, 
                            metrics=model_metrics)
        model.fit(X_train, 
                  y_train, 
                  epochs=model_epochs, 
                  batch_size=model_batch_size, 
                  class_weight=class_weight,
                  validation_split=0.2,
                  verbose=1)

        # evaluate the model
        print("Evaluate...")
        result = model.evaluate(X_test, y_test, verbose=0)
        result_dict = dict(zip(model.metrics_names, result))
        # Store the performance metrics for each repetition, store both train and test results
        for key, value in result_dict.items():
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
del results
tf.keras.backend.clear_session()
