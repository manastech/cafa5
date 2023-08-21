import warnings
warnings.filterwarnings("ignore")

from manas_cafa.models.singleTermConv1DModel import singleTermConv1DModel
from manas_cafa.bio.protein import Protein

import pandas as pd
from Bio import SeqIO
import random
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt
import argparse
import os

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Some constants
GO_TERMS_TRAIN_SET = '../data/go_terms_train_set_maxlen500_minmembers75.tsv'
CANDIDATE_TERMS_TEST_SET = '../data/go_terms_test_candidates_maxlen500_minmembers75.tsv'
TRAIN_SEQUENCES_PATH = '../cafa-5-protein-function-prediction/Train/train_sequences.fasta'
UNIPROT_ENTRIES_PATH = '/data/uniprot/entries/'

def singleTermPipeline( go_term, epochs, batch_size, save_path = ""):
    print("RUNNING PIPELINE FOR GO term:", go_term)

    # Read the train set
    train_set_df = pd.read_csv(GO_TERMS_TRAIN_SET, index_col=0, sep='\t', header=None, names=['go_term', 'proteins'] )
    candidate_row = train_set_df.loc[go_term]

    # Get the proteins for the go term
    candidate_row_proteins = candidate_row.proteins.split(',')

    # Load the proteins (max length 500)
    other_proteins = []
    for record in SeqIO.parse(TRAIN_SEQUENCES_PATH, 'fasta'):
        if len(record.seq) <= 500:
            other_proteins.append(record.id)
    other_proteins = list(set(other_proteins) - set(candidate_row_proteins))

    # select random proteins from the other proteins
    random.shuffle(other_proteins)
    other_proteins = other_proteins[:min(len(candidate_row_proteins)*10,15000)]

    # get the one-hot encoded sequence of the proteins in the candidate row
    candidate_row_proteins_one_hot = []
    for protein_id in candidate_row_proteins:
        protein = Protein(protein_id)
        protein.load_file(UNIPROT_ENTRIES_PATH+protein_id[-2:]+"/"+protein_id+".xml")
        candidate_row_proteins_one_hot.append(protein.one_hot_sequence())
    
    # fix the length of the one-hot encoded sequences to 500
    candidate_row_proteins_one_hot = [np.pad(seq, ((0, 500-len(seq)), (0,0)), 'constant') for seq in candidate_row_proteins_one_hot if len(seq) <= 500]

    # get the one-hot encoded sequence of the other proteins
    other_proteins_one_hot = []
    for protein_id in other_proteins:
        protein = Protein(protein_id)
        try:
            protein.load_file(UNIPROT_ENTRIES_PATH+protein_id[-2:]+"/"+protein_id+".xml")
            other_proteins_one_hot.append(protein.one_hot_sequence())
        except:
            continue

    # fix the length of the one-hot encoded sequences to 500
    other_proteins_one_hot = [np.pad(seq, ((0, 500-len(seq)), (0,0)), 'constant') for seq in other_proteins_one_hot if len(seq) <= 500]

    print("Positive examples:", len(candidate_row_proteins_one_hot))
    print("Negative examples:", len(other_proteins_one_hot))

    # create singleTermConv1DModel object
    model = singleTermConv1DModel(input_shape=(500, 20), num_classes=2)

    # compile the model
    model.compile_model()

    # create the train and test sets
    eighty_percent_candidate = int(len(candidate_row_proteins_one_hot)*0.9)
    eighty_percent_other = int(len(other_proteins_one_hot)*0.9)

    x_train = np.array(candidate_row_proteins_one_hot[:eighty_percent_candidate] + other_proteins_one_hot[:eighty_percent_other])
    y_train = np.array([[1,0]]*eighty_percent_candidate + [[0,1]]*eighty_percent_other)
    x_test = np.array(candidate_row_proteins_one_hot[eighty_percent_candidate:] + other_proteins_one_hot[eighty_percent_other:])
    y_test = np.array([[1,0]]*(len(candidate_row_proteins_one_hot)-eighty_percent_candidate) + [[0,1]]*(len(other_proteins_one_hot)-eighty_percent_other))

    # fit the model, silently
    model.fit_model(x_train, y_train, x_test, y_test,
                epochs=epochs, batch_size=batch_size, verbose=0)
    
    # save the model to a file in /data/models/<go_term[-2:]>/<go_term>, create the directory if it doesn't exist
    if save_path != "":
        model.model.save(save_path)

    # plot a roc curve for the test set
    y_pred_keras = model.predict(x_test).ravel()
    fpr_keras, tpr_keras, thresholds_keras = roc_curve(y_test.ravel(), y_pred_keras)
    auc_keras = auc(fpr_keras, tpr_keras)

    print("AUC:", auc_keras)
    print("")

if __name__ == "__main__":
    # ask for terminal parameters and call the singleTermPipeline function
    parser = argparse.ArgumentParser(description='Train a model for a single GO term.')
    parser.add_argument('go_term', metavar='go_term', type=str, help='The GO term to train the model for.')
    parser.add_argument('--epochs', metavar='epochs', type=int, default=10, help='The number of epochs to train the model for.')
    parser.add_argument('--batch_size', metavar='batch_size', type=int, default=32, help='The batch size to train the model for.')
    parser.add_argument('--save_path', metavar='save_path', type=int, default=32, help='The path to save the model to (optional)')
    args = parser.parse_args()
    singleTermPipeline(args.go_term, args.epochs, args.batch_size, args.save_path)
