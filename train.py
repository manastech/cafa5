import tensorflow as tf
import pandas as pd
import numpy as np
from tensorflow.keras import layers, models
from Bio import SeqIO
import csv
from collections import defaultdict

TRAIN_SEQUENCES_PATH = "./cafa-5-protein-function-prediction/Train/train_sequences.fasta"
TARGET_TERMS_PATH = './MFO_level_2_terms.tsv'

ALPHABET = list('ACDEFGHIKLMNPQRSTVWY')
MAX_PROTEIN_LENGTH = 500 # In this toy model we will consider proteins up to 500 amino acids long

def go_terms_and_proteins(terms_path, protein_list):
    proteins_vs_terms = defaultdict(lambda:[])
    go_terms = []
    with open(terms_path) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        for protein, go_term in rd:
            if protein in protein_list:
                proteins_vs_terms[protein].append(go_term)
                go_terms.append(go_term)

    go_terms = list(set(go_terms))
    return go_terms, proteins_vs_terms
 
def encode_seq_onehot(sequence):
    char_to_int = dict((c, i) for i, c in enumerate(ALPHABET))
    integer_encoded = [char_to_int[char] if char in ALPHABET else 0 for char in sequence]
    onehot_encoded = list()

    for value in integer_encoded:
        letter = [0 for _ in range(len(ALPHABET))]
        letter[value] = 1
        onehot_encoded.append(letter)

    return np.array(onehot_encoded) 

def encode_go_terms_in_proteins(go_terms, go_terms_per_protein, proteins_ids):
    encoded_go_terms_per_protein = [np.zeros(len(go_terms)) for _ in range(len(proteins_ids))]

    for i, protein in enumerate(proteins_ids):
        for go_term in go_terms_per_protein[protein]:
            encoded_go_terms_per_protein[i][go_terms.index(go_term)] = 1
            
    return np.array(encoded_go_terms_per_protein)

with open(TRAIN_SEQUENCES_PATH) as handle:
    training_set_ids = [record.id for record in SeqIO.parse(handle, "fasta") if len(record.seq) < MAX_PROTEIN_LENGTH ]

with open(TRAIN_SEQUENCES_PATH) as handle:
    encoded_training_seqs = np.array([
        encode_seq_onehot((record.seq+"X"*MAX_PROTEIN_LENGTH)[0:MAX_PROTEIN_LENGTH]) 
        for record in SeqIO.parse(handle, "fasta") 
        if len(record.seq) < MAX_PROTEIN_LENGTH
    ])

print("Total proteins", len(training_set_ids))
print("Encoded input shape:", encoded_training_seqs.shape)

go_terms_mfo_level_2, proteins_with_mfo_level_2_terms = go_terms_and_proteins(TARGET_TERMS_PATH, training_set_ids)

print('GO terms of interest: ', len(go_terms_mfo_level_2))
print('Proteins having a GO term of interest: ', len(proteins_with_mfo_level_2_terms))

encoded_go_terms = encode_go_terms_in_proteins(go_terms_mfo_level_2, proteins_with_mfo_level_2_terms, training_set_ids)

print("Encoded output shape:", encoded_go_terms.shape)

model = models.Sequential()
model.add(layers.Conv1D(64, (64), activation='relu', input_shape=(500,20)))
model.add(layers.MaxPooling1D())
model.add(layers.Conv1D(32, (32), activation='relu'))
model.add(layers.MaxPooling1D())
model.add(layers.Flatten())
model.add(layers.Dense(163, activation = 'sigmoid'))

model.compile(
    loss = 'binary_crossentropy', 
    optimizer = "adam",
    metrics = ['accuracy']
)

print(model.summary())
model.fit(encoded_training_seqs, encoded_go_terms, batch_size= 50,epochs=10)

model.save('./OneHotEncodingSeqModel.h5')

