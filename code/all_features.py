'''
Created on Dec 7, 2015

@author: Nicolas Metts

This module runs classification on feature sets for the large yeast data set.

'''
import argparse
import itertools
import classify

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

ALL_CLASSIFIERS = ['log_reg', 'svm', 'ada_boost', 'knn', 'gnb']

def __classify_all_features(main_args):
    """
    Runs classification or kfold validation on the specified feature sets

    Args:
        main_args: Arguments to specify train/test files, classification, and
                    other arguments relative to the task
    """
    feature_sets = []
    for i in range(args.features_min - 1, main_args.features_max, 1):
        combinations = itertools.combinations(ALL_FEATURES, i + 1)
        for combination in combinations:
            feature_sets.append(list(combination))
    print "Length of all combinations: " + str(len(feature_sets))
    for feature_set in feature_sets:
        classify_args = classify.ClassifyArgs(features=feature_set,
                                               classifiers=ALL_CLASSIFIERS,
                                               write_to_log=True,
                                               train_file=main_args.train_file,
                                               test_file=main_args.test_file,
                                               data_file=main_args.data_file,
                                               classify=main_args.classify,
                                               kfold=main_args.kfold)
        classify.main(classify_args)
        classify_args = classify.ClassifyArgs(features=feature_set,
                                               classifiers=ALL_CLASSIFIERS,
                                               write_to_log=True,
                                               train_file=main_args.train_file,
                                               test_file=main_args.test_file,
                                               data_file=main_args.data_file,
                                               vote=main_args.vote,
                                               classify=main_args.classify,
                                               kfold=main_args.kfold)
        classify.main(classify_args)
    print "feature set loop done"

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--classify", help="Classify using training and test set",
                           action="store_true")
    argparser.add_argument("--kfold", help="10-fold cross validation",
                           action="store_true")
    argparser.add_argument("--feature_min", help="Minimum number of features",
                           type=int, default=49)
    argparser.add_argument("--feature_max", help="Maximum number of features",
                           type=int, default=50)
    argparser.add_argument("--data_file", help="Name of data file",
                           type=str, default="../data/large_yeast_data.csv",
                           required=False)
    argparser.add_argument("--train_file", help="Name of train file",
                           type=str, default="../data/testing_data_large.csv",
                           required=False)
    argparser.add_argument("--test_file", help="Name of test file",
                           type=str, default="../data/training_data_large.csv",
                           required=False)
    argparser.add_argument("--vote",
                           help="Ensemble classifier. 'hard' = majority, 'soft' = average",
                           type=str, default='none')
    args = argparser.parse_args()
    __classify_all_features(args)

