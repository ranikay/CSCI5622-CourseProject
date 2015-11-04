'''
Created on Nov 3, 2015

@author: Nicolas Metts

This module is responsible for selecting features to be used in classification
and for classifying and measuring the accuracy of classifications
'''

import argparse
from csv import DictReader
import numpy as np

ORF = 'ORF'
MITOCHONDRIA = 'mitochondria'
CYTOPLASM = 'cytoplasm'
ER = 'er'
NUCLEUS = 'nucleus'
VACUOLE = 'vacuole'
OTHER = 'other'
CAI = 'CAI'
NC = 'NC'
GC = 'GC'
L_GA = 'L_aa'
GRAVY = 'Gravy'
DOV_EXPR = 'DovEXPR'
BLAST_HITS_IN_YEAST = 'BLAST_hits_in_yeast'
INTXN_PARTNERS = 'intxn_partners'
CHROMOSOME = 'chromosome'
CHR_POSITION = 'chr_position'
INTRON = 'intron'
CLOSE_STOP_RATIO = 'close_stop_ratio'
RARE_AA_RATIO = 'rare_aa_ratio'
TM_HELIX = 'tm_helix'
IN_HOW_MANY_OF_5_PROKS = 'in_how_many_of_5_proks'
IN_HOW_MANY_OF_6_CLOSE_YEAST = 'in_how_many_of_6_close_yeast'
SGD_ESS = 'SGD_ess'

def create_feature_file(source, features, out_file_name, train_file):
    """
    A function to create a feature file

    Args:
        source: The data source
        features: An iterable collection of features
        out_file_name: The name of the output file
        train_file: A boolean indicating whether this is a training file or not
    """
    # TODO: Implement this!
    data = []
    for example in source:
        row = []
        if not train_file:
            row.append(example[SGD_ESS])
        for feature in features:
            # For now, just appending the features. May need to do some 
            # pre-processing on the raw data
            row.append(example[feature])
        data.append(row)
    return np.array(data)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--train_file", help="Name of train file",
                           type=str, default="../data/training_data.csv", required=False)
    argparser.add_argument("--test_file", help="Name of test file",
                           type=str, default="../data/testing_data.csv", required=False)
    argparser.add_argument("--create_features", help="Create a training file",
                           action="store_true")
    argparser.add_argument("--classify", help="Classify using training and test set",
                           action="store_true")
    argparser.add_argument("--classifiers", help="A list of classifiers to use",
                           nargs='+', required=False)
    argparser.add_argument("--metrics", help="A list of metrics to use",
                           nargs='+', required=False)
    # Is this option needed if we're using training and test files?
    argparser.add_argument("--cross_validate", help="Cross validate using training and test set",
                           action="store_true")
    argparser.add_argument("--features", help="Features to be used",
                           nargs='+', required=False)
    args = argparser.parse_args()
    # Cast to list to keep it all in memory
    train = list(DictReader(open(args.train_file, 'r')))
    test = list(DictReader(open(args.test_file, 'r')))