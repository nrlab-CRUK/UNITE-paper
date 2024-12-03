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
from sklearn.metrics import confusion_matrix

# make wd and x_file as parameter to this script using argparse

parser = argparse.ArgumentParser(description="Train CNN model ",
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

parser.add_argument("--model_file", 
                    help="model file name with absolute path", 
                    required=True, 
                    type=str, 
                    dest="model_file")





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
performance_file = 'performance' + args.x_file[7:-4] + '.binary_model.csv'
performance_file_path = os.path.join(args.wd, performance_file)

# load the data
x = np.load(x_file_path, allow_pickle=True)
# remove the 188th row from 2nd dimension
x = np.delete(x, 276, axis = 1)

# Remove the first 50 rows from the 3rd dimension of x
#x = np.delete(x, np.s_[0:50], axis = 2)


y = np.load(y_file_path, allow_pickle=True)
z = np.load(z_file_path, allow_pickle=True)

y = np.array(y)


# print tensor shape
print('x shape: ', x.shape)
print('y shape: ', y.shape)
print('z shape: ', z.shape)


# print the count of each class
unique, counts = np.unique(y, return_counts=True)
print('unique: ', unique)
print('counts: ', counts)


# print model file full path

print('model_file: ', args.model_file)

# metrics



# define model
def create_model_simpler(input_shape, loss, learning_rate, metrics):
    model = keras.Sequential(
        [
            keras.Input(shape=input_shape),

            layers.Conv2D(8, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(3, 2)),

            layers.Conv2D(8, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),

            layers.Conv2D(16, kernel_size=(3, 3), activation="relu"),
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

################################################################################
# load the existing model
################################################################################
print("Loading model...")
# load the model and custom metrics
model = keras.models.load_model(args.model_file, custom_objects={'sens': sens, 'spec': spec, 'get_f1': get_f1})

# evaluate the model and save the results


# evaluate the model using test data
print("Evaluating...")
test_result = model.evaluate(x, y, verbose=2)
# add "test_" to the beginning of each metric
test_metrics = ['test_' + s for s in model.metrics_names]
# create dictionary of metric names and results
test_result_dict = dict(zip(test_metrics, test_result))

print(test_result_dict)


# create a pandas DataFrame from the new dictionary, one column contaiing the metric names and the other the results
df = pd.DataFrame(list(test_result_dict.items()), columns=['metric', 'result'])
#df = pd.DataFrame.from_dict(test_result_dict, orient='index', columns=)

# save to csv
df.to_csv(performance_file_path, index=False)
print('Results saved to: ', performance_file_path)
print(df)
print('Done!')


# get the confusion matrix and print it
y_pred = model.predict(x)
y_pred = np.round(y_pred)
cm = confusion_matrix(y, y_pred)
# extract the true negatives, false positives, false negatives, and true positives
tn, fp, fn, tp = cm.ravel()
# print the confusion matrix, true negative, false positive, false negative, and true positive
print('Confusion Matrix:')
print(cm)
print('True Negative: ', tn)
print('False Positive: ', fp)
print('False Negative: ', fn)
print('True Positive: ', tp)

