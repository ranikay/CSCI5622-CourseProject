'''
Created on Dec 7, 2015

@author: Nicolas Metts
'''
import classify
import sys
from numpy.random import choice
import itertools

ALL_FEATURES = ["Transcript.length", "Strand", "GC.content", "Enzyme", 
                "SEG.low.complexity", "Transmembrane.domain", "Signal.peptide",
                "Coiled.coil", "Nucleus", "Mitochondria", "ER", "Cytoplasm",
                "Ribosome", "Interaction.partners", "Mw", "PI", "Protein.length",
                "Gravy.score", "Aromaticity.Score", "CAI", "Codon.Bias",
                "FOP.Score", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His",
                "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser",
                "Thr", "Val", "Trp", "Tyr", "Carbons", "Hydrogens", "Nitrogens",
                "Oxygens", "Sulphurs", "Instability.index", "Alphatic.index",
                "EC.number"]

ALL_CLASSIFIERS = ['log_reg', 'svm', 'ada_boost', 'knn','gnb']

if __name__ == '__main__':
    features_min = int(sys.argv[1])
    features_max = int(sys.argv[2])
    feature_sets = []
    """
    for i in range(features_min, features_max + 1, 1):
        print "Getting combinations for " + str(i) + " features"
        chosen = choice(ALL_FEATURES, i, replace=False)
        for c in chosen:
            feature_sets.append(list(c))"""
    for i in range(features_min - 1, features_max, 1):
        chosen = itertools.combinations(ALL_FEATURES, i + 1)
        for c in chosen:
            feature_sets.append(list(c))
    print "Length of all combinations: " + str(len(feature_sets))
    for feature_set in feature_sets:
        args = classify.classify_args(features=feature_set,
                                      classifiers=ALL_CLASSIFIERS,
                                      write_to_log=True, 
                                      train_file="../data/testing_data_large.csv",
                                      test_file="../data/training_data_large.csv",
                                      data_file='large_yeast_data.csv',
                                      classify=True)
        classify.main(args)
        args = classify.classify_args(features=feature_set,
                                      classifiers=ALL_CLASSIFIERS,
                                      write_to_log=True, 
                                      train_file="../data/testing_data_large.csv",
                                      test_file="../data/training_data_large.csv",
                                      data_file='large_yeast_data.csv',
                                      vote='hard',
                                      classify=True)
        classify.main(args)
