'''
Created on Nov 3, 2015

@author: Nicolas Metts

This module is responsible for selecting features to be used in classification
and for classifying and measuring the accuracy of classifications
'''

import argparse
from csv import DictReader
import numpy as np
import random
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import SGDClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.dummy import DummyClassifier
from sklearn.svm import SVC
from sklearn import cross_validation
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import precision_score, recall_score, accuracy_score
import matplotlib.pyplot as plt
from itertools import compress

# Constants for feature names
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
#added features
CHRM_AND_POS = 'chrm_and_pos' #chromosome and postion 
LOC_AND_5_PROKS = 'loc_and_5_proks' #combine location and how many of 5 proks
YEAST_AND_PROKS = 'yeast_and_proks' #combine the number of matches in yeast and proks
#label
SGD_ESS = 'SGD_ess'


#Constants for classifier names
LOG_REG = 'log_reg'
SVM = 'svm'
ADA_BOOST = 'ada_boost'
KNN = 'knn'
GNB = 'gnb' #Gaussian Naive Bayes
UNIFORM = 'uniform' #DummyClassifier
# TODO: Add more classifiers

def write_log(args, classifier, accuracy, precision, recall,
              true_count, actual_count, X_train, X_test):
    """
    Function to write results of a run to a file.
    """

    # Get the kernel type if classifier is SVM, otherwise just put NA
    get_kernel = lambda x: x == 'svm' and args.kernel or "NA"

    # Log important info
    with open('../logs/log_table.txt', 'a') as f: f.write("\n")
    log = [str(args.data_file),str(args.train_file),str(args.test_file),str(args.create_features),
           str(classifier), str(get_kernel(classifier)),str(args.scale), str(len(X_train)), str(len(X_test)),
           str(precision), str(recall), str(accuracy), str(true_count), str(actual_count), str([args.features])]
    line = "\t".join(log)

    with open('../logs/log_table.txt', 'a') as f:
        f.write(line)

    
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

def pairwise_graphs(data):
    """
    This function produces graphs for each pair of features given a list of dictionaries
    A LOT of graphs are produced. I recommend constraining the keys in the loops
    if using show()

    Args:
        data: a list of dictionaries
    """
    #first remove ORF
    for e in data:
        e.pop(ORF, None)

    keys = data[0].keys()
    labels = [int(e[SGD_ESS]) for e in data]
    notLabels = [0 if l else 1 for l in labels]
    i = 0
    for x in keys:
        X = [float(e[x]) for e in data]
        xPos = list(compress(X, labels))
        xNeg = list(compress(X, notLabels))
        i += 1
        for y in keys[i:]:
            Y = [float(e[y]) for e in data]
            yPos = list(compress(Y, labels))
            yNeg = list(compress(Y, notLabels))

            pos = plt.scatter(xPos, yPos, c='r', alpha=0.5)
            neg = plt.scatter(xNeg, yNeg, c='g', alpha=0.5)
            title = x+' vs '+y
            plt.title(title)
            plt.xlabel(x)
            plt.ylabel(y)
            plt.legend((pos, neg), ('Essential', 'Non-Essential'), scatterpoints=1)
            #plt.show()
            plt.savefig('../data/graphs/'+title+'.png')

def single_var_graphs(data):
    """
    This function is for visualizing single variables against labels

    Args:
        data: a list of dictionaries
    """
    #first remove ORF
    for e in data:
        e.pop(ORF, None)

    keys = data[0].keys()
    labels = [int(e[SGD_ESS]) for e in data]
    notLabels = [0 if l else 1 for l in labels]
    i = 0
    for x in keys[:5]: #
        X = [float(e[x]) for e in data]
        xPos = list(compress(X, labels))
        xNeg = list(compress(X, notLabels))
        
        title = x+' Essentiality'
        axes = plt.gca()
        axes.set_ylim([-1,2])
        pos = plt.scatter(xPos, [1] * len(xPos), c='r', alpha=0.5)
        neg = plt.scatter(xNeg, [0] * len(xNeg), c='g', alpha=0.5)
        plt.title(title)
        plt.xlabel(x)
        plt.ylabel('Essentiality')
        plt.legend((pos, neg), ('Essential', 'Non-Essential'), scatterpoints=1)
        plt.show()
        #plt.savefig('../data/graphs/'+title+'.png')

def svm_classify(train_X, train_Y, test_X, test_Y, kernel, reg):
    """
    A function to run an SVM classification

    Args:
        train_X: training feature values
        train_Y: training labels
        test_X: testing feature values
        test_Y: testing labels
        kernel: a string representing the kernel to use
        reg: a float representing the regularization parameter
    """

    clf = SVC(kernel=kernel, C=reg)
    clf.fit(train_X, train_Y)
    sc = clf.score(test_X, test_Y)
    print('SVM score', kernel, reg, sc)

    return clf

def __print_and_log_results(clf, x_test, y_test):
    predictions = clf.predict(x_test)
    accuracy = accuracy_score(y_test, predictions, [0, 1])
    precision = precision_score(y_test, predictions, [0, 1])
    recall = recall_score(y_test, predictions, [0, 1])
    print "Train/test set sizes: " + str(len(X_train)) + "/" + str(len(x_test))
    print "Precision is: " + str(precision)
    print "Recall is: " + str(recall)
    print "Accuracy is: " + str(accuracy)
    true_count = len([1 for p in predictions if p=='1'])
    actual_count = len([1 for y in y_test if y=='1'])
    print "True count (prediction/actual): " + str(true_count) + "/" + str(actual_count)
    if args.write_to_log:
    # Write out results as a table to log file
        write_log(out_file_name = "../logs/log_table.txt", args = args,
                    classifier = classifier, accuracy = accuracy,
                    precision = precision, recall = recall,
                    true_count = true_count, actual_count = actual_count,
                    X_train = X_train, X_test = X_test)

def __get_classifier_model(classifier, args):
    """
    Convenience function for obtaining a classification model

    Args:
        classifier (str): A string indicating the name of the classifier
        args: An arguments object

    Returns:
        A classification model based on the given classifier string
    """
    # Make SGD Logistic Regression model the default
    model = SGDClassifier(loss='log', penalty='l2', shuffle=True)
    if classifier == LOG_REG:
        model = SGDClassifier(loss='log', penalty='l2', shuffle=True)
    elif classifier == SVM:
        model = SVC(kernel=args.kernel)
    elif classifier == ADA_BOOST:
        model = AdaBoostClassifier()
    elif classifier == KNN:
        model = KNeighborsClassifier(n_neighbors=5, algorithm='ball_tree')
    elif classifier == GNB:
        model = GaussianNB()
    elif classifier == UNIFORM:
        model = DummyClassifier()
    return model

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--data_file", help="Name of data file",
                           type=str, default="../data/cerevisiae_compiled_features.csv", 
                           required=False)
    argparser.add_argument("--train_file", help="Name of train file",
                           type=str, default="../data/training_data.csv", required=False)
    argparser.add_argument("--test_file", help="Name of test file",
                           type=str, default="../data/testing_data.csv", required=False)
    argparser.add_argument("--create_features", help="Create a training file",
                           action="store_true")
    argparser.add_argument("--classify", help="Classify using training and test set",
                           action="store_true")
    argparser.add_argument("--classifiers", help="A list of classifiers to use",
                           nargs='+', required=False, default = ['svm'])
    argparser.add_argument("--metrics", help="A list of metrics to use",
                           nargs='+', required=False)
    argparser.add_argument("--kernel", help="The kernel to be used for SVM classification",
                           type=str, default='rbf')
    # Is this option needed if we're using training and test files?
    argparser.add_argument("--cross_validate", help="Cross validate using training and test set",
                           action="store_true")
    argparser.add_argument("--write_to_log", help="Send output to log file",
                           action="store_true")
    argparser.add_argument("--features", help="Features to be used",
                           nargs='+', required=True)
    argparser.add_argument("--scale", help="Scale the data with StandardScale",
                           action="store_true")
    args = argparser.parse_args()
    

    if args.classify:
        # Cast to list to keep it all in memory
        train = list(DictReader(open(args.train_file, 'r')))
        test = list(DictReader(open(args.test_file, 'r')))

        labels = []
        for line in train:
            if not line[SGD_ESS] in labels:
                labels.append(line[SGD_ESS])

        train_features = []
        for example in train:
            train_feat = []
            for feature in args.features:
                train_feat.append(example[feature])
            train_features.append(train_feat)
        x_train = np.array(train_features, dtype=float)
        
        test_features = []
        for example in test:
            test_feature = []
            for feature in args.features:
                test_feature.append(example[feature])
            test_features.append(test_feature)
        x_test = np.array(test_features, dtype=float)
        
        y_train = np.array(list(labels.index(x[SGD_ESS])
                         for x in train))
        
        y_test = np.array(list(labels.index(x[SGD_ESS])
                         for x in test))

        if args.scale:
            scaler = StandardScaler()
            x_train = scaler.fit_transform(x_train)
            x_test = scaler.fit_transform(x_test)

        # Train classifier
        for classifier in args.classifiers:
            # TODO: Add some way to specify options for classifiers
            model = __get_classifier_model(classifier, args)
            clf = model.fit(x_train, y_train)
            print "Using classifier " + classifier
            __print_and_log_results(clf, x_test, y_test)
    
    elif args.cross_validate:
        # Cast to list to keep it all in memory
        data = list(DictReader(open(args.data_file, 'rU')))
        labels = []
        for line in data:
            labels.append(line[SGD_ESS])
        train_features = []
        for example in data:
            train_feat = []
            for feature in args.features:
                train_feat.append(example[feature])
            train_features.append(train_feat)
        x_train = np.array(train_features, dtype=float)
        rand_int = random.randint(1, len(data))
        X_train, X_test, y_train, y_test = cross_validation.train_test_split \
            (x_train,labels, test_size=0.1, random_state=13)
        if args.scale:
            scaler = StandardScaler()
            X_train = scaler.fit_transform(X_train)
            X_test = scaler.fit_transform(X_test)
        for classifier in args.classifiers:
            model = __get_classifier_model(classifier, args)
            clf = model.fit(X_train, y_train)
            print "Using classifier " + classifier
            __print_and_log_results(clf, X_test, y_test)

