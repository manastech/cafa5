import warnings
warnings.filterwarnings("ignore")

from manas_cafa.bio.protein import Protein

import pandas as pd
from Bio import SeqIO
import numpy as np
from keras.models import load_model
import tensorflow as tf
import tensorflow.keras as keras
import argparse
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Some constants
TEST_SEQUENCES = '../cafa-5-protein-function-prediction/Test/test_sequences.fasta'
UNIPROT_ENTRIES_PATH = '/data/uniprot/entries/'


def singleTermPredPipeline(go_term, results_path, batch_size=2500):
    print("PREDICTING", go_term)
    # load the model
    model = load_model(go_term)

    # Load the proteins (max length 500)
    proteins = []
    for record in SeqIO.parse(TEST_SEQUENCES, 'fasta'):
        if len(record.seq) <= 500:
            proteins.append(record.id)

    # predict the proteins in batches
    for i in range(0, len(proteins), batch_size):
        batch = proteins[i:i+batch_size]
        proteins_one_hot = []
        
        for protein_id in batch:
            try:
                protein = Protein(protein_id)
                protein.load_file(UNIPROT_ENTRIES_PATH+protein_id[-2:]+"/"+protein_id+".xml")
                proteins_one_hot.append(protein.one_hot_sequence())
            except:
                pass

        print("Predicting for", len(proteins_one_hot), "proteins")

        # fix the length of the one-hot encoded sequences to 500
        proteins_one_hot = [np.pad(seq, ((0, 500-len(seq)), (0,0)), 'constant') for seq in proteins_one_hot if len(seq) <= 500]

        # predict the probability of the proteins in the candidate row
        predictions = model.predict(np.array(proteins_one_hot))

        # get the probability of the proteins in the candidate row
        proteins_prob = [p[0] for p in predictions]

        # get an ordered list of the proteins in the candidate row with their probabilities
        proteins_prob_ordered = list(zip(proteins, proteins_prob))
        proteins_prob_ordered.sort(key=lambda x: x[1], reverse=True)

        # save the predictions to a file in CAFA5 format
        with open(results_path, 'a') as f:
            for p in proteins_prob_ordered:
                #round the probability to 3 decimals
                if round(p[1], 3) > 0.0:
                    f.write(p[0]+'\tGO:'+go_term.name.replace(".h5", "")+'\t'+str(round(p[1], 3))+'\n')

if __name__ == "__main__":
    # ask for terminal parameters and call the singleTermPipeline function
    parser = argparse.ArgumentParser(description='Train a model for a single GO term.')
    parser.add_argument('go_term', metavar='go_term', type=str, help='The GO term to train the model for.')
    args = parser.parse_args()
    singleTermPredPipeline(args.go_term)
