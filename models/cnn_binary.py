import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np
from sklearn.utils import class_weight
import keras.backend as K
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
import os
import argparse
import datetime
from tensorflow.keras.callbacks import TensorBoard
# make wd and x_file as parameter to this script using argparse

parser = argparse.ArgumentParser(description="Cross Validation for CNN model ",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--wd", 
                    help="working directory", 
                    type=str, 
                    dest="wd", 
                    default='.')

parser.add_argument("--x_file", 
                    help="x file name", 
                    required=True, 
                    type=str, 
                    dest="x_file")

parser.add_argument("--num_repeats", 
                    help="n repeats for CV", 
                    type=int, 
                    dest="num_repeats",
                    default=20)

parser.add_argument("--num_folds", 
                    help="n folds for CV", 
                    type=int, 
                    dest="num_folds",
                    default=5)

parser.add_argument("--model_epochs",
                    help="number of epochs",
                    type=int,
                    dest="model_epochs",
                    default=35)

parser.add_argument("--model_batch_size",
                    help="batch size",
                    type=int,
                    dest="model_batch_size",
                    default=16)

parser.add_argument("--model_verbose",
                    help="verbose",
                    type=int,
                    dest="model_verbose",
                    default=1)

parser.add_argument("--model_validation_split",
                    help="validation split",
                    type=float,
                    dest="model_validation_split",
                    default=0.3)

parser.add_argument("--model_learning_rate",
                    help="learning rate",
                    type=float,
                    dest="model_learning_rate",
                    default=0.0008)




args = parser.parse_args()

#set working dir
x_file_path = os.path.join(args.wd, args.x_file)
y_file = 'y' + args.x_file[1:]
y_file_path = os.path.join(args.wd, y_file)
z_file = 'z' + args.x_file[1:]
z_file_path = os.path.join(args.wd, z_file)

# print file names
print('x_file_path: ', x_file_path)

# set output file names
performance_csv = 'performance' + args.x_file[7:-4] + '.debug.csv'
performance_csv_path = os.path.join(args.wd, performance_csv)

# load the data
x = np.load(x_file_path, allow_pickle=True)
y = np.load(y_file_path, allow_pickle=True)
z = np.load(z_file_path, allow_pickle=True)
y = np.array(y)

# print tensor shape
print('x shape: ', x.shape)
print('y shape: ', y.shape)
print('z shape: ', z.shape)

# compute class weights 
model_class_weight = dict(zip(np.unique(y), class_weight.compute_class_weight(class_weight='balanced',
                                                                        classes=np.unique(y),
                                                                        y=y)))


# tensor board callback
log_folder = "logs/fit/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
tensorboard_callback = TensorBoard(log_dir=log_folder,
                         histogram_freq=1,
                         write_graph=True,
                         write_images=True,
                         update_freq='epoch',
                         profile_batch=2,
                         embeddings_freq=1)

model_callbacks = [tensorboard_callback]


# cv parameters
num_repeats = args.num_repeats
num_folds = args.num_folds
kfold_grouped = StratifiedGroupKFold(n_splits=num_folds,
                                     shuffle=True)

# model parameters
model_input_shape = x.shape[1:]
model_learning_rate = 0.0008
model_loss = 'binary_crossentropy'

# training parameters
model_epochs = args.model_epochs
model_batch_size = args.model_batch_size 
model_verbose = args.model_verbose 
model_validation_split = args.model_validation_split 


# define model
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
            # layers.Dropout(0.2),
            layers.Dense(1, activation="sigmoid")
        ])

    opt = keras.optimizers.Adam(learning_rate=learning_rate)
    model.compile(loss=loss,
                  optimizer=opt,
                  metrics=metrics)
    return model


# define function to calculate metrics

def sens(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    return true_positives / (possible_positives + K.epsilon())


def spec(y_true, y_pred):
    true_negatives = K.sum(K.round(K.clip((1 - y_true) * (1 - y_pred), 0, 1)))
    possible_negatives = K.sum(K.round(K.clip(1 - y_true, 0, 1)))
    return true_negatives / (possible_negatives + K.epsilon())

def get_f1(y_true, y_pred): #taken from old keras source code
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    recall = true_positives / (possible_positives + K.epsilon())
    f1_val = 2*(precision*recall)/(precision+recall+K.epsilon())
    return f1_val


model_metrics = [
    'accuracy',
    tf.keras.metrics.AUC(curve="ROC", name='AUROC'),
    tf.keras.metrics.AUC(curve="PR", name='AUPRC'),
    tf.keras.metrics.Precision(name='PPV'),
    tf.keras.metrics.Recall(name='recall_sens'),
    # add f1 score
    get_f1,
    #tf.keras.metrics.F1Score(name='f1'),
    tf.keras.metrics.SensitivityAtSpecificity(specificity=0.95, name="sen_95spe"),
    tf.keras.metrics.SensitivityAtSpecificity(specificity=0.98, name="sen_98spe"),
    tf.keras.metrics.SensitivityAtSpecificity(specificity=0.99, name="sen_99spe")
]

###############################################################################
# perform CV
###############################################################################

results = {}

for i in range(num_repeats):
    # report which loop we are on
    print('# ' * 50)
    print('Repeat: ', i + 1)

    # split into train and test
    for train_idx, test_idx in kfold_grouped.split(x, y, groups=z):
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
                  class_weight=model_class_weight,
                  validation_split=model_validation_split,
                  verbose=model_verbose,
                  callbacks=model_callbacks)

        # evaluate the model using test data
        print("Evaluate...")
        test_result = model.evaluate(X_test, y_test, verbose=0)
        # add "test_" to the beginning of each metric
        test_metrics = ['test_' + s for s in model.metrics_names]
        # create dictionary of metric names and results
        test_result_dict = dict(zip(test_metrics, test_result))

        # evaluate the model using train data
        train_result = model.evaluate(X_train, y_train, verbose=0)
        # add "train_" to the beginning of each metric
        train_metrics = ['train_' + s for s in model.metrics_names]
        # create dictionary of metric names and results
        train_result_dict = dict(zip(train_metrics, train_result))

        # concatenate the train and test dictionaries
        result_dict = {**train_result_dict, **test_result_dict}


        # Store the performance metrics for each repetition 
        for key, value in result_dict.items():
            if key not in results:
                results[key] = []
            results[key].append(value)

# Create a pandas DataFrame from the results dictionary

# remove 'estimator' key
relevant_keys = [k for k in results.keys() if k != 'estimator']
filtered_dict = {k: v for k, v in results.items() if k in relevant_keys}

# create a pandas DataFrame from the new dictionary
cv_df = pd.DataFrame.from_dict(filtered_dict)
# save to csv
cv_df.to_csv(performance_csv_path, index=False)
print('Results saved to: ', performance_csv_path)
print(cv_df)
print('Done!')

# clean up memory
del results
tf.keras.backend.clear_session()