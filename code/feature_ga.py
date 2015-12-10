# use a genetic algorithm to find good combinations of features with respect to some ML technique

import argparse
from csv import DictReader
import numpy as np
import random
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics.classification import precision_score, recall_score, accuracy_score
from sklearn.preprocessing import StandardScaler
from itertools import compress
import matplotlib.pyplot as plt

# Constants for feature names
ORF = 'ORF'
MITOCHONDRIA = 'Mitochondria'
CYTOPLASM = 'Cytoplasm'
ER = 'ER'
NUCLEUS = 'Nucleus'
VACUOLE = 'Vacuole'
OTHER = 'Other.location'
CAI = 'CAI'
NC = 'NC'
GC = 'GC.content'
L_GA = 'Protein.length'
GRAVY = 'Gravy.score'
DOV_EXPR = 'DovEXPR'
BLAST_HITS_IN_YEAST = 'BLAST.hits_in_yeast'
INTXN_PARTNERS = 'Interaction.partners'
CHROMOSOME = 'Chromosome'
CHR_POSITION = 'Chr.position'
INTRON = 'Intron'
CLOSE_STOP_RATIO = 'Close.stop.ratio'
RARE_AA_RATIO = 'Rare.aa_ratio'
TM_HELIX = 'TM.helix'
IN_HOW_MANY_OF_5_PROKS = 'in_how_many_of_5_proks'
IN_HOW_MANY_OF_6_CLOSE_YEAST = 'in_how_many_of_6_close_yeast'
# added features
CHRM_AND_POS = 'chrm_and_pos' # chromosome and postion
LOC_AND_5_PROKS = 'loc_and_5_proks' # combine location and how many of 5 proks
YEAST_AND_PROKS = 'yeast_and_proks' # combine the number of matches in yeast and proks
# label
ESSENTIAL = 'Essential'


def init_pop(n, num_features):
    return [np.random.randint(0, 2, num_features) for _ in xrange(n)]

def calc_fit(model, metric, train_x, train_y, test_x, test_y, p):
    train_x = map(lambda x: list(compress(x, p)), train_x)
    test_x = map(lambda x: list(compress(x, p)), test_x)
    clf = model.fit(train_x, train_y)
    predictions = clf.predict(test_x)
    if metric == 'precision': return precision_score(test_y, predictions, [0, 1])
    elif metric == 'recall': return recall_score(test_y, predictions, [0, 1])
    elif metric == 'accuracy': return accuracy_score(test_y, predictions, [0, 1])
    return precision_score(test_y, predictions, [0, 1]) + recall_score(test_y, predictions, [0, 1]) + accuracy_score(test_y, predictions, [0, 1])

# simple single point crossover. we could try other types
def recombine(pop, fits):
    # include the top two individual with no changes
    new_pop = list(zip(* sorted(zip(fits, pop), key=lambda x: x[0], reverse=True)[:2])[1])
    while len(new_pop) < len(pop):
        # pick individuals
        (x, y) = np.random.choice(range(len(pop)), size=2, p=fits)
        a = pop[x]
        b = pop[y]
        # pick an index
        i = np.random.randint(0, len(pop[0]))
        # add new individuals
        new_pop.insert(0, np.concatenate([a[:i], b[i:]]))
        new_pop.insert(0, np.concatenate([b[:i], a[i:]]))
    return new_pop

def mutate(pop, rate):
    # pick an individual to mutate and then an index to flip
    new_pop = []
    for p in pop:
        if rate > np.random.random():
            i = np.random.randint(0, len(p))
            d = (p[i] + 1) % 2 # 0 if p[i] == 1 else 0
            np.put(p, i, d)
            new_pop.insert(0, p)
        else: new_pop.insert(0, p)
    return new_pop

#
def feature_graphs(all_feature_scores):
    small_data = {}
    large_data = {}

    # this is fairly specific to the exact format of these files. they could be made more general if necessary
    with open('../logs/small_feature_ga_run_log.txt', 'r') as small:
        for m in ['precision', 'accuracy', 'recall', 'all']:
            small_data[m] = []
            small.readline() # info line
            while True:
                l = small.readline().strip()
                print('line', l)
                if l == '***': break
                small_data[m].append(float(l))
            small.readline() # empty line
    with open('../logs/large_feature_ga_run_log.txt', 'r') as large:
        for m in ['precision', 'accuracy', 'recall', 'all']:
            large_data[m] = []
            large.readline() # info line
            while True:
                l = large.readline().strip()
                if l == '***': break
                large_data[m].append(float(l))
            large.readline() # empty line

    # graphs
    # small
    plt.plot(range(100), small_data['accuracy'], 'r-', label='Accuracy, all: ' + str(all_feature_scores['small']['accuracy'])[:5] + ', max: ' + str(max(small_data['accuracy']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['accuracy']] * 100, 'r--') # , label='Accuracy all features')
    plt.plot(range(100), small_data['precision'], 'b-', label='Precision, all: ' + str(all_feature_scores['small']['precision'])[:5] + ', max: ' + str(max(small_data['precision']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['precision']] * 100, 'b--') # , label='Precision all features')
    plt.plot(range(100), small_data['recall'], 'g-', label='Recall, all: ' + str(all_feature_scores['small']['recall'])[:5] + ', max: ' + str(max(small_data['recall']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['recall']] * 100, 'g--') # , label='Recall all features')
    axes = plt.gca()
    axes.set_ylim([0, 1.])
    title = 'Small Data Set: Fitness Over Time'
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(loc=4)
    plt.savefig('../data/graphs/' + title + '.png')
    plt.show()

    # large
    plt.plot(range(100), large_data['accuracy'], 'r-', label='Accuracy, all: ' + str(all_feature_scores['large']['accuracy'])[:5] + ', max: ' + str(max(large_data['accuracy']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['accuracy']] * 100, 'r--') # , label='Accuracy all features')
    plt.plot(range(100), large_data['precision'], 'b-', label='Precision, all: ' + str(all_feature_scores['large']['precision'])[:5] + ', max: ' + str(max(large_data['precision']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['precision']] * 100, 'b--') # , label='Precision all features')
    plt.plot(range(100), large_data['recall'], 'g-', label='Recall, all: ' + str(all_feature_scores['large']['recall'])[:5] + ', max: ' + str(max(large_data['recall']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['recall']] * 100, 'g--') # , label='Recall all features')
    axes = plt.gca()
    axes.set_ylim([0, 1.])
    title = 'Large Data Set: Fitness Over Time'
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(loc=4)
    plt.savefig('../data/graphs/' + title + '.png')
    plt.show()

    # Small vs large comparisons
    # accuracy
    plt.plot(range(100), large_data['accuracy'], 'b-', label='Accuracy Large, all: ' + str(all_feature_scores['large']['accuracy'])[:5] + ', max: ' + str(max(large_data['accuracy']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['accuracy']] * 100, 'b--') # , label='Accuracy all features')
    plt.plot(range(100), small_data['accuracy'], 'r-', label='Accuracy Small, all: ' + str(all_feature_scores['small']['accuracy'])[:5] + ', max: ' + str(max(small_data['accuracy']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['accuracy']] * 100, 'r--') # , label='Accuracy all features')
    axes = plt.gca()
    axes.set_ylim([.7, .85])
    title = 'Small vs Large Data Set: Accuracy Over Time'
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(loc=4)
    plt.savefig('../data/graphs/' + title + '.png')
    plt.show()

    # precision
    plt.plot(range(100), large_data['precision'], 'b-', label='Precision Large, all: ' + str(all_feature_scores['large']['precision'])[:5] + ', max: ' + str(max(large_data['precision']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['precision']] * 100, 'b--') # , label='Precision all features')
    plt.plot(range(100), small_data['precision'], 'r-', label='Precision, all: ' + str(all_feature_scores['small']['precision'])[:5] + ', max: ' + str(max(small_data['precision']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['precision']] * 100, 'r--') # , label='Precision all features')
    axes = plt.gca()
    axes.set_ylim([.35, .6])
    title = 'Small vs Large Data Set: Precision Over Time'
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(loc=4)
    plt.savefig('../data/graphs/' + title + '.png')
    plt.show()

    # recall
    plt.plot(range(100), large_data['recall'], 'b-', label='Recall Large, all: ' + str(all_feature_scores['large']['recall'])[:5] + ', max: ' + str(max(large_data['recall']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['recall']] * 100, 'b--') # , label='Recall all features')
    plt.plot(range(100), small_data['recall'], 'r-', label='Recall Small, all: ' + str(all_feature_scores['small']['recall'])[:5] + ', max: ' + str(max(small_data['recall']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['recall']] * 100, 'r--') # , label='Recall all features')
    axes = plt.gca()
    axes.set_ylim([0, .4])
    title = 'Small vs Large Data Set: Recall Over Time'
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(loc=4)
    plt.savefig('../data/graphs/' + title + '.png')
    plt.show()

    # all
    plt.plot(range(100), large_data['all'], 'b-', label='All Large, all: ' + str(all_feature_scores['large']['all'])[:5] + ', max: ' + str(max(large_data['all']))[:5])
    plt.plot(range(100), [all_feature_scores['large']['all']] * 100, 'b--') # , label='Recall all features')
    plt.plot(range(100), small_data['all'], 'r-', label='All Small, all: ' + str(all_feature_scores['small']['all'])[:5] + ', max: ' + str(max(small_data['all']))[:5])
    plt.plot(range(100), [all_feature_scores['small']['all']] * 100, 'r--') # , label='Recall all features')
    axes = plt.gca()
    axes.set_ylim([1.3, 1.65])
    title = 'Small vs Large Data Set: All Over Time'
    plt.title(title)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.legend(loc=4)
    plt.savefig('../data/graphs/' + title + '.png')
    plt.show()

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--classifier', help='the ml technique to use',
                                                 type=str, default='knn')
    argparser.add_argument('--metric', help='fitness metric. precision, recall, accuracy, or all',
                                                 type=str, default='precision')
    argparser.add_argument("--data_file", help="Name of data file",
                                                 type=str, default="../data/large_yeast_data.csv", required=False)
    argparser.add_argument("--train_file", help="Name of train file",
                                                 type=str, default="../data/testing_data_large.csv", required=False) # flip train and test. should be ~ 80/20
    argparser.add_argument("--test_file", help="Name of test file",
                                                 type=str, default="../data/training_data_large.csv", required=False)
    argparser.add_argument('--scale', help="Scale the data with StandardScale",
                                                 action="store_true")
    argparser.add_argument('--m_rate', help='Mutation rate',
                                                 type=float, default=0.01)
    argparser.add_argument('--generations', help='Number of generations',
                                                 type=int, default=100)
    argparser.add_argument('--pop_size', help='Population size',
                                                 type=int, default=20)
    argparser.add_argument('--write_gens', help='Append information from this run to ../log/feature_ga_run_log.txt',
                                                 action="store_true")
    args = argparser.parse_args()

    print('Read in data')
    train = list(DictReader(open(args.train_file, 'r')))
    test = list(DictReader(open(args.test_file, 'r')))

    # generate additional features to add and test
    # these can be added to the data file later
    # combine chromosome and location: CHROMOSOME + CHR_POSITION
    # chrm_lengths = [230218, 813184, 316620, 1531933, 576874, 270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333, 1091291, 948066, 85779, 6318] #chrm 1-16, mito, 2-micron.. how are the last two reprented numerically? are they relevant at all?
    # location and how many of 5 proks: mitochondria = 10, cytoplasm = 20, er = 30, nucelus = 40, vacuoule = 50, otehr = 60
    # loc_dic = {MITOCHONDRIA: 10, CYTOPLASM: 20, ER: 30, NUCLEUS: 40, VACUOLE: 50, OTHER: 60}
    # how many hits in yeast and proks

    # add new features to test and train
    # currently only in small data set: chrm, chrm position, vacuole, other, how many proks
    features = []
    if args.train_file == '../data/small_yeast_data.csv':
        for d in train:
            d[CHRM_AND_POS] = float(d[CHROMOSOME]) + float(d[CHR_POSITION])
            d[LOC_AND_5_PROKS] = list(compress([10, 20, 30, 40, 50, 60], [d[e] for e in [MITOCHONDRIA, CYTOPLASM, ER, NUCLEUS, VACUOLE, OTHER]]))[0] + float(d[IN_HOW_MANY_OF_5_PROKS])
            d[YEAST_AND_PROKS] = float(d[IN_HOW_MANY_OF_6_CLOSE_YEAST]) + float(d[IN_HOW_MANY_OF_5_PROKS])
        for d in test:
            d[CHRM_AND_POS] = float(d[CHROMOSOME]) + float(d[CHR_POSITION])
            d[LOC_AND_5_PROKS] = list(compress([10, 20, 30, 40, 50, 60], [d[e] for e in [MITOCHONDRIA, CYTOPLASM, ER, NUCLEUS, VACUOLE, OTHER]]))[0] + float(d[IN_HOW_MANY_OF_5_PROKS])
            d[YEAST_AND_PROKS] = float(d[IN_HOW_MANY_OF_6_CLOSE_YEAST]) + float(d[IN_HOW_MANY_OF_5_PROKS])
        # ORF not included
        features = [MITOCHONDRIA, CYTOPLASM, ER, NUCLEUS, VACUOLE, OTHER, CAI, NC, GC, L_GA, GRAVY, DOV_EXPR, BLAST_HITS_IN_YEAST, INTXN_PARTNERS, CHROMOSOME, CHR_POSITION, INTRON, CLOSE_STOP_RATIO, RARE_AA_RATIO, TM_HELIX, IN_HOW_MANY_OF_5_PROKS, IN_HOW_MANY_OF_6_CLOSE_YEAST, CHRM_AND_POS, LOC_AND_5_PROKS, YEAST_AND_PROKS]
    else:
        # remove NA entries completely
        train = filter(lambda d: 'NA' not in d.values(), train)
        test = filter(lambda d: 'NA' not in d.values(), test)
        # also remove ORF if present
        features = [f for f in train[0].keys() if f != 'ORF']

    for d in test:
        if d[ESSENTIAL] != '1' and d[ESSENTIAL] != '0': print(d, len(d.keys()))
    train_x = [[float(e[f]) for f in e if f != ESSENTIAL and f != ORF] for e in train]
    train_y = [float(e[ESSENTIAL]) for e in train]
    test_x = [[float(e[f]) for f in e if f != ESSENTIAL and f != ORF] for e in test]
    test_y = [float(e[ESSENTIAL]) for e in test]

    #
    if args.scale:
        print('Scale data')
        scaler = StandardScaler()
        train_x = scaler.fit_transform(train_x)
        test_x = scaler.fit_transform(test_x)

    print('Set model')
    model = None
    if args.classifier == 'knn': model = KNeighborsClassifier(n_neighbors=5, algorithm='ball_tree')

    gen_maxes = []

    # initialize the population
    print('Initialize Population')
    pop = init_pop(args.pop_size, len(features))

    print('evolve')
    m = (0, pop[0])
    for g in xrange(args.generations):
        print('Generation', g)
        # determine fitnesses
        fits = [calc_fit(model, args.metric, train_x, train_y, test_x, test_y, p) for p in pop]
        fitTot = sum(fits)
        normedFits = [float(x) / fitTot for x in fits]
        print('max fit', max(fits))
        gen_maxes.append(max(fits))
        # update max
        mIndex = fits.index(max(fits))
        m = (max(fits), pop[mIndex]) if max(fits) > m[0] else m
        # Recombination
        pop = recombine(pop, normedFits)
        # Mutation
        pop = mutate(pop, args.m_rate)

    print('Result', args.metric)
    pop_with_fits = sorted(zip([calc_fit(model, args.metric, train_x, train_y, test_x, test_y, p) for p in pop], pop), key=lambda x: x[0], reverse=True)
    for (f, p) in pop_with_fits:
        print(p, f)
    print('max encountered', m[1], list(compress(features, m[1])), m[0])
    if args.metric == 'all':
        print('max precision: ', calc_fit(model, 'precision', train_x, train_y, test_x, test_y, m[1]))
        print('max accuracy: ', calc_fit(model, 'accuracy', train_x, train_y, test_x, test_y, m[1]))
        print('max recall: ', calc_fit(model, 'recall', train_x, train_y, test_x, test_y, m[1]))
    print('all features', 'precision', calc_fit(model, 'precision', train_x, train_y, test_x, test_y, [1] * len(features)))
    print('all features', 'accuracy', calc_fit(model, 'accuracy', train_x, train_y, test_x, test_y, [1] * len(features)))
    print('all features', 'recall', calc_fit(model, 'recall', train_x, train_y, test_x, test_y, [1] * len(features)))
    print('all features', 'all', calc_fit(model, 'all', train_x, train_y, test_x, test_y, [1] * len(features)))
    print('Done')

    all_feature_scores = {'small':{}, 'large':{}}
    all_feature_scores['small']['precision'] = 0.51449275362318836
    all_feature_scores['small']['accuracy'] = 0.78421409214092141
    all_feature_scores['small']['recall'] = 0.22015503875968992
    all_feature_scores['small']['all'] = 1.5188618845237998
    all_feature_scores['large']['precision'] = 0.44327176781002636
    all_feature_scores['large']['accuracy'] = 0.80085261875761271
    all_feature_scores['large']['recall'] = 0.17910447761194029
    all_feature_scores['large']['all'] = 1.4232288641795794
    # feature_graphs(all_feature_scores)

    if args.write_gens:
        with open('../logs/large_feature_ga_run_log.txt', 'a') as F:
            F.write(args.metric + ',' + str(args.generations) + ',' + str(args.pop_size) + ',' + str(args.m_rate) + '\n') # run info
            for e in gen_maxes:
                F.write(str(e) + '\n')
            F.write('***\n')
