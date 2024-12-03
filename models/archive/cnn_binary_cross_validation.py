
# see https://github.com/tirthajyoti/Deep-learning-with-Python/blob/master/Notebooks/Keras_Scikit_Learn_wrapper.ipynb
#see https://towardsdatascience.com/are-you-using-the-scikit-learn-wrapper-in-your-keras-deep-learning-model-a3005696ff38

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

from sklearn.preprocessing import minmax_scale

import numpy as np
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['figure.dpi']=150



def create_model(input_shape,
metrics=['accuracy'], 
loss='binary_crossentropy', 
learning_rate=0.001):
    model = keras.Sequential(
        [
            keras.Input(shape=input_shape),
            
            layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(3, 2)),
            
            layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),
            
            layers.Conv2D(32, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),
            
            layers.Conv2D(16, kernel_size=(3, 3), activation="relu"),
            layers.MaxPooling2D(pool_size=(2, 2)),
            
            layers.Flatten(),
            layers.Dense(32, activation="relu"),
            #layers.Dropout(0.2),
            layers.Dense(1, activation="sigmoid")
        ])
    
    opt = keras.optimizers.Adam(learning_rate=learning_rate)
    model.compile(loss=loss,
                    optimizer=opt,
                    metrics= metrics)
    return model



  

from tensorflow.keras.wrappers.scikit_learn import KerasClassifier

model = KerasClassifier(build_fn=create_model, 
                        epochs=10, 
                        batch_size=32, 
                        verbose=0)
                        

from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score

num_folds = 10
kfold = StratifiedKFold(n_splits=num_folds, 
                        shuffle=True)

cv_results = cross_val_score(model, 
                          x, y, 
                          cv=kfold, scoring= metrics,
                          verbose=2)
                          
