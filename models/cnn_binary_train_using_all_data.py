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
from sklearn.model_selection import train_test_split

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

parser.add_argument("--model_type", 
                    help="simpler or deeper", 
                    required=True, 
                    type=str, 
                    default='simpler',
                    dest="model_type")

parser.add_argument("--model_epochs",
                    help="number of epochs",
                    type=int,
                    dest="model_epochs",
                    default=35)

parser.add_argument("--model_batch_size",
                    help="batch size",
                    type=int,
                    dest="model_batch_size",
                    default=32)

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
model_type = args.model_type

# model file name with parameters and timestamp
model_file = 'model' + args.x_file[7:-4] + '.binary_model.' + model_type + '.epochs_' + str(args.model_epochs) + '.batch_size_' + str(args.model_batch_size) + '.verbose_' + str(args.model_verbose) + '.validation_split_' + str(args.model_validation_split) + '.learning_rate_' + str(args.model_learning_rate) + '.h5'
# add path to model file
model_file_path = os.path.join(args.wd, model_file)


# load the data
x = np.load(x_file_path, allow_pickle=True)
x = np.delete(x, 276, axis = 1)
#x = np.delete(x, np.s_[0:50], axis = 2)


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

# print class weights
print('class weights: ', model_class_weight)


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






# Create and train the model
if args.model_type == 'simpler':
     model = create_model_simpler(input_shape=model_input_shape,
                     loss=model_loss,
                     learning_rate=model_learning_rate,
                     metrics=model_metrics)
else:
     model = create_model(input_shape=model_input_shape,
                     loss=model_loss,
                     learning_rate=model_learning_rate,
                     metrics=model_metrics)


# split data into train and validation, keep class balance
x_train, x_val, y_train, y_val = train_test_split(x, y, test_size=model_validation_split, stratify=y)


# print horizontal line
print('----------------------------------------')

# print data shape
print('x_train shape:', x_train.shape)
print('x_val shape:', x_val.shape)
print('y_train shape:', y_train.shape)
print('y_val shape:', y_val.shape)

# print class balance
print('Class balance in train set:', np.unique(y_train, return_counts=True))
print('Class balance in validation set:', np.unique(y_val, return_counts=True))

# print horizontal line
print('----------------------------------------')

# print model summary
print(model.summary())

# print horizontal line
print('----------------------------------------')


model.fit(x_train,
          y_train,
          epochs=model_epochs,
          batch_size=model_batch_size,
          class_weight=model_class_weight,
          validation_data=(x_val, y_val),
          verbose=model_verbose,
          callbacks=model_callbacks)


# save model
model.save(model_file_path)
print('Model saved to:', model_file_path)

