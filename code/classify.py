'''
Created on Nov 3, 2015

@author: Nicolas Metts

This module is responsible for selecting features to be used in classification
and for classifying and measuring the accuracy of classifications
'''

import argparse
from csv import DictReader
import csv
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
from sklearn.ensemble import VotingClassifier
#from mlxtend.classifier import EnsembleClassifier
from sklearn.metrics import precision_score, recall_score, accuracy_score
import matplotlib.pyplot as plt
from itertools import compress

# Constants for feature names

ORF = 'ORF'
#label
ESSENTIAL = 'Essential'


#Constants for classifier names
LOG_REG = 'log_reg'
SVM = 'svm'
ADA_BOOST = 'ada_boost'
KNN = 'knn'
GNB = 'gnb' #Gaussian Naive Bayes
UNIFORM = 'uniform' #DummyClassifier
# TODO: Add more classifiers

class classify_args:
    def __init__(self, data_file="../data/small_yeast_data.csv",
                 train_file="../data/training_data.csv", 
                 test_file="../data/testing_data.csv",create_features=False,
                 classify=False, classifiers=['svm'], kernel = 'rbf',
                 cross_validate=False, write_to_log=False, features=[],
                 scale=False, vote='none'):
        self.data_file = data_file
        self.train_file = train_file
        self.test_file = test_file
        self.create_features = create_features
        self.classify = classify
        self.classifiers = classifiers
        self.kernel = kernel
        self.cross_validate = cross_validate
        self.write_to_log = write_to_log
        self.features = features
        self.scale = scale
        self.vote = vote

def write_log(out_file_name, args, classifier, accuracy, precision, recall,
              true_count, actual_count, X_train, X_test, all_features):
    """
    Function to write results of a run to a file.
    """

    # Get the kernel type if classifier is SVM, otherwise just put NA
    get_kernel = lambda x: x == 'svm' and args.kernel or "NA"

    # Log important info
    log = [args.data_file,args.train_file,args.test_file,
           classifier, get_kernel(classifier),args.scale, len(X_train),
           len(X_test), precision, recall, accuracy, true_count, actual_count]
    # Include a TRUE/FALSE column for each feature
    log += [feature in args.features for feature in all_features]

    with open(out_file_name, 'a') as f:
        out_writer = csv.writer(f, lineterminator='\n')
        out_writer.writerow(log)

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
            row.append(example[ESSENTIAL])
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
    labels = [int(e[ESSENTIAL]) for e in data]
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
    labels = [int(e[ESSENTIAL]) for e in data]
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

def __print_and_log_results(clf, classifier, x_train, x_test, y_test, out_file_name,
                            all_features, args):
    predictions = clf.predict(x_test)
    accuracy = accuracy_score(y_test, predictions, [0, 1])
    precision = precision_score(y_test, predictions, [0, 1])
    recall = recall_score(y_test, predictions, [0, 1])
    print "Train/test set sizes: " + str(len(x_train)) + "/" + str(len(x_test))
    print "Precision is: " + str(precision)
    print "Recall is: " + str(recall)
    print "Accuracy is: " + str(accuracy)
    true_count = len([1 for p in predictions if p==1])
    actual_count = len([1 for y in y_test if y==1])
    print "True count (prediction/actual): " + str(true_count) + "/" + str(actual_count)
    if args.write_to_log:
    # Write out results as a table to log file
        write_log(out_file_name=out_file_name, args=args, classifier=classifier,
                    accuracy=accuracy, precision=precision, recall=recall,
                    true_count=true_count, actual_count=actual_count,
                    X_train=x_train, X_test=x_test, all_features=all_features)

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
    if(args.vote == 'none'):
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
    else:
        #We might consider passing all individual classifiers back to compare to the ensemble
        #See the last line in http://scikit-learn.org/stable/modules/ensemble.html#id24
        clfs = []
        for clf in args.classifiers:
            if clf == LOG_REG:
                clfs.append((clf,SGDClassifier(loss='log', penalty='l2', shuffle=True)))
            elif clf == SVM:
                clfs.append((clf,SVC(kernel=args.kernel)))
            elif clf == ADA_BOOST:
                clfs.append((clf, AdaBoostClassifier()))
            elif clf == KNN:
                clfs.append((clf,KNeighborsClassifier(n_neighbors=5, algorithm='ball_tree')))
            elif clf == GNB:
                clfs.append((clf, GaussianNB()))
            elif clf == UNIFORM:
                clfs.append((clf, DummyClassifier()))
        model = VotingClassifier(estimators=clfs, voting=args.vote)

    return model


def main(args):
    voting_methods = ['none', 'hard', 'soft']
    assert args.vote in voting_methods, "--vote must be one of 'none', 'hard', 'soft'"
    
    if 'small_yeast_data' in args.data_file:
        out_file_name = '../logs/small_yeast_data_log.txt'
    if 'large_yeast_data' in args.data_file:
        out_file_name = '../logs/large_yeast_data_log.txt'

    if args.classify:
        # Store column names as features, except ORF and Essential
        all_features = DictReader(open(args.train_file, 'r')).fieldnames
        all_features.remove('ORF')
        all_features.remove('Essential')
        
        # Cast to list to keep it all in memory
        train = list(DictReader(open(args.train_file, 'r')))
        test = list(DictReader(open(args.test_file, 'r')))

        labels = []
        for line in train:
            if not line[ESSENTIAL] in labels:
                labels.append(line[ESSENTIAL])

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
        
        y_train = np.array(list(labels.index(x[ESSENTIAL])
                         for x in train))
        
        y_test = np.array(list(labels.index(x[ESSENTIAL])
                         for x in test))

        if args.scale:
            scaler = StandardScaler()
            x_train = scaler.fit_transform(x_train)
            x_test = scaler.fit_transform(x_test)
        
        if args.vote == 'none':
            for classifier in args.classifiers:
                model = __get_classifier_model(classifier, args)
                clf = model.fit(x_train, y_train)
                print "Using classifier " + classifier
                __print_and_log_results(clf, classifier, x_train, x_test, y_test,
                                        out_file_name, all_features, args)
        else:
            model = __get_classifier_model('none',args)
            clf = model.fit(x_train, y_train)
            print "Using classifier: vote "+ args.vote + " with ", args.classifiers
            classifier = "vote-" + args.vote + "-with-classifiers_"
            classifier += "_".join(args.classifiers)
            __print_and_log_results(clf, classifier, x_train, x_test, y_test,
                                    out_file_name, all_features, args)
    
    elif args.cross_validate:

        # Store column names as features, except ORF and Essential
        all_features = DictReader(open(args.data_file, 'rU')).fieldnames
        all_features.remove('ORF')
        all_features.remove('Essential')
        # Cast to list to keep it all in memory
        data = list(DictReader(open(args.data_file, 'rU')))
        
        labels = []
        for line in data:
            labels.append(int(line[ESSENTIAL]))
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
            scaler = StandardScaler().fit(X_train)
            X_train = scaler.transform(X_train)
            X_test = scaler.transform(X_test)
        if args.vote == 'none':
            for classifier in args.classifiers:
                model = __get_classifier_model(classifier, args)
                clf = model.fit(X_train, y_train)
                print "Using classifier " + classifier
                __print_and_log_results(clf, classifier, X_train, X_test, y_test,
                                        out_file_name, all_features, args)
        else:
            model = __get_classifier_model('none',args)
            clf = model.fit(X_train, y_train)
            print "Using classifier: vote "+ args.vote + " with ", args.classifiers
            classifier = "vote-" + args.vote + "-with-classifiers_"
            classifier += "_".join(args.classifiers)
            __print_and_log_results(clf, classifier, X_train, X_test, y_test, out_file_name,
                                    all_features, args)
    
if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--data_file", help="Name of data file",
                           type=str, default="../data/small_yeast_data.csv", 
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
    argparser.add_argument("--vote", help="Ensemble classifier. 'hard' = majority, 'soft' = average",
                           type=str, default='none')
    args = argparser.parse_args()
    
    main(args)
