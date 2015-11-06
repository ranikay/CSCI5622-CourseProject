#use a genetic algorithm to find good combinations of features with respect to some ML technique

import argparse
from csv import DictReader
import numpy as np
import random
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics.classification import precision_score, recall_score, accuracy_score
from sklearn.preprocessing import StandardScaler
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
SGD_ESS = 'SGD_ess'


def init_pop(n, num_features):
  return [np.random.randint(0,2,num_features) for _ in xrange(n)]

def calc_fit(model, metric, train_x, train_y, test_x, test_y, p):
  train_x = map(lambda x: list(compress(x, p)), train_x)
  test_x = map(lambda x: list(compress(x, p)), test_x)
  clf = model.fit(train_x, train_y)
  predictions = clf.predict(test_x)
  if metric == 'precision': return precision_score(test_y, predictions, [0,1])
  elif metric == 'recall': return recall_score(test_y, predictions, [0,1])
  return accuracy_score(test_y, predictions, [0,1])

#simple single point crossover. we could try other types
def recombine(pop, fits):
  #include the top two individual with no changes
  new_pop = list(zip(* sorted(zip(fits, pop), key=lambda x: x[0], reverse=True)[:2])[1])
  while len(new_pop) < len(pop):
    #pick individuals
    (x, y) = np.random.choice(range(len(pop)), size=2, p=fits)
    a = pop[x]
    b = pop[y]
    #pick an index
    i = np.random.randint(0,len(pop[0]))
    #add new individuals
    new_pop.insert(0, np.concatenate([a[:i], b[i:]]))
    new_pop.insert(0, np.concatenate([b[:i], a[i:]]))
  return new_pop

def mutate(pop, rate):
  #pick an individual to mutate and then an index to flip
  new_pop = []
  for p in pop:
    if rate > np.random.random():
      i = np.random.randint(0, len(p))
      d = (p[i]+1) % 2 #0 if p[i] == 1 else 0 
      np.put(p, i, d)
      new_pop.insert(0, p)
    else: new_pop.insert(0, p)
  return new_pop
      

if __name__ == '__main__':
  argparser = argparse.ArgumentParser()
  argparser.add_argument('--classifier', help='the ml technique to use',
                         type=str, default='knn')
  argparser.add_argument('--metric', help='fitness metric. precision, recall, or accuracy',
                         type=str, default='precision')
  argparser.add_argument("--data_file", help="Name of data file",
                         type=str, default="../data/cerevisiae_compiled_features.csv", required=False)
  argparser.add_argument("--train_file", help="Name of train file",
                         type=str, default="../data/training_data.csv", required=False)
  argparser.add_argument("--test_file", help="Name of test file",
                         type=str, default="../data/testing_data.csv", required=False)
  argparser.add_argument('--scale', help="Scale the data with StandardScale",
                         action="store_true")
  argparser.add_argument('--m_rate', help='Mutation rate',
                         type=float, default=0.01)
  argparser.add_argument('--generations', help='Number of generations',
                         type=int, default=100)
  argparser.add_argument('--pop_size', help='Population size',
                         type=int, default=20)
  args = argparser.parse_args()


  print('Read in data')
  train = list(DictReader(open(args.train_file, 'r')))
  test = list(DictReader(open(args.test_file, 'r')))

  train_x = [[float(e[f]) for f in e if f != SGD_ESS and f != ORF] for e in train]
  train_y = [float(e[SGD_ESS]) for e in train]
  test_x = [[float(e[f]) for f in e if f != SGD_ESS and f != ORF] for e in test]
  test_y = [float(e[SGD_ESS]) for e in test]

  #ORF not included
  features = [MITOCHONDRIA, CYTOPLASM, ER, NUCLEUS, VACUOLE, OTHER, CAI, NC, GC, L_GA, GRAVY, DOV_EXPR, BLAST_HITS_IN_YEAST, INTXN_PARTNERS, CHROMOSOME, CHR_POSITION, INTRON, CLOSE_STOP_RATIO, RARE_AA_RATIO, TM_HELIX, IN_HOW_MANY_OF_5_PROKS, IN_HOW_MANY_OF_6_CLOSE_YEAST]

  if args.scale:
    print('Scale data')
    scaler = StandardScaler()
    train_x = scaler.fit_transform(train_x)
    test_x = scaler.fit_transform(test_x)

  print('Set model')
  model = None
  if args.classifier == 'knn': model = KNeighborsClassifier(n_neighbors=5, algorithm='ball_tree')

  #initialize the population
  print('Initialize Population')
  pop = init_pop(args.pop_size, len(features))

  print('evolve')
  m = (0, pop[0])
  for g in xrange(args.generations):
    print('Generation', g)
    #determine fitnesses
    fits = [calc_fit(model, args.metric, train_x, train_y, test_x, test_y, p) for p in pop]
    fitTot = sum(fits)
    normedFits = [float(x)/fitTot for x in fits]
    print('max fit', max(fits))
    #update max
    mIndex = fits.index(max(fits))
    m = (max(fits), pop[mIndex]) if max(fits) > m[0] else m
    #Recombination
    pop = recombine(pop, normedFits)
    #Mutation
    pop = mutate(pop, args.m_rate)

  print('Result', args.metric)
  pop_with_fits = sorted(zip([calc_fit(model, args.metric, train_x, train_y, test_x, test_y, p) for p in pop], pop), key=lambda x: x[0], reverse=True)
  for (f, p) in pop_with_fits:
    print(p, f)
  print('max encountered', m[1], m[0])
  print('all features', 'precision', calc_fit(model, 'precision', train_x, train_y, test_x, test_y, [1]*len(features)))
  print('all features', 'accuracy', calc_fit(model, 'accuracy', train_x, train_y, test_x, test_y, [1]*len(features)))
  print('all features', 'recall', calc_fit(model, 'recall', train_x, train_y, test_x, test_y, [1]*len(features)))
  print('Done')





