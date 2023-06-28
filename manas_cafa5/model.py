import tensorflow as tf
import numpy as np
from tensorflow.keras.preprocessing.sequence import pad_sequences
from random import sample

class Model:
    def __init__(self, model, max_protein_length):
        self.model = model
        self.max_protein_length = max_protein_length

    @staticmethod
    def cnn1d(layers, max_protein_length=500):
        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Conv1D(32, kernel_size=3, activation='relu', input_shape=(max_protein_length,20)))
        for layer in layers:
            model.add(layer)
        model.add(tf.keras.layers.Dense(2, activation='softmax')) # output layer
        return Model(model, max_protein_length)

    # todo: subclass tf.keras.Model
    def call(self, inputs):
        return self.model(inputs)

    @staticmethod
    def parse_trainingset(file):
        terms = {}
        with open(file, 'r', encoding='utf-8') as f:
            for line in f:
                [go_term,associated_uniprots] = line.rstrip().split('\t')
                terms[go_term] = associated_uniprots.split(',')
        return terms

    def compile(self, trainingset, graph, protein_uniprot_id, protein_cache, max_distance=1):
        protein = protein_cache.load(protein_uniprot_id)
        children = protein.go_terms_children(graph, max_distance)

        uniprot_ids = trainingset[protein_uniprot_id]
        proteins = { uid: protein_cache.load(uid) for uid in uniprot_ids }

        x_train = [
            pad_sequences(
                protein.one_hot_sequence(),
                maxlen=self.max_protein_length,
                padding='post',
                truncating='post'
            )
            for uid in uniprot_ids
        ]
        y_train = [ uid in children for uid in uniprot_ids ]
        return LoadedModel(
            model=self.model,
            inputs=x_train,
            outputs=y_train,
        )

class LoadedModel:
    def __init__(self, model, inputs, outputs):
        self.model = model
        self.inputs = inputs
        self.outputs = outputs

    def fit(self, **kwargs):
        kwargs = kwargs.copy()
        kwargs.update({ 'inputs': self.inputs, 'outputs': self.outputs })
        return self.model.fit(**kwargs)

#def rnn1d(x_train, y_train, x_val, y_val, layers, max_protein_length=500):
#    pass
