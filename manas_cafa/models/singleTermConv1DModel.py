# A class to model in Keras for single term prediction using Conv1D
# taking as input a set of proteins one hot encoded

import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

class singleTermConv1DModel( object ):
    
    def __init__(self, input_shape, num_classes):
        self.input_shape = input_shape
        self.num_classes = num_classes
        self.model = self.create_model()

    def create_model(self):
        # Create the model
        model = keras.Sequential()
        model.add(layers.Conv1D(32, kernel_size=3, activation='relu', input_shape=self.input_shape))
        model.add(layers.Conv1D(64, kernel_size=3, activation='relu'))
        model.add(layers.Dropout(0.5))
        model.add(layers.MaxPooling1D(pool_size=2))
        model.add(layers.Flatten())
        model.add(layers.Dense(128, activation='relu'))
        model.add(layers.Dense(self.num_classes, activation='softmax'))
        return model

    def compile_model(self):
        self.model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['accuracy'])

    def fit_model(self, x_train, y_train, x_val, y_val, epochs, batch_size, verbose=1):
        self.model.fit(x_train, y_train, validation_data=(x_val, y_val), epochs=epochs, batch_size=batch_size, verbose=verbose)

    def evaluate_model(self, x_test, y_test):
        self.model.evaluate(x_test, y_test)

    def predict(self, x):
        return self.model.predict(x)




