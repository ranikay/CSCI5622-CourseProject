#Train our set of classifiers on some subset of features using the yeast data in large data testing file (flipped training and testing sizes)
#then use these classifier on human data

import classify
from csv import DictReader
import numpy as np
from itertools import compress

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

FEATURE_SUBSETS = [[True,True,True,True,True,True,True,False,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True,True]]

ALL_CLASSIFIERS = ['log_reg', 'svm', 'ada_boost', 'knn','gnb']

#focus on mcf7 as the label gene... essential when 2 or more stds below avg
#MCF7_BREAST is index 116
def gen_human_data_file(features, labels, out):
  human_features_a = list(DictReader(open(features, 'r'), delimiter='\t'))
  human_labels = list(DictReader(open(labels, 'r'), delimiter='\t'))

  #remove features that are not in the list of all features
  human_features = []
  for entry in human_features_a:
    keys = entry.keys()
    d = {}
    for k in keys:
      #added a few keys that I believe are synonymous
      #replace names with those used in yeast data
      if k == 'SEG': d['SEG.low.complexity'] = entry[k]
      if k == 'Avg.GC.content': d['GC.content'] = entry[k]
      if k == 'Mitochondrion': d['Mitochondria'] = entry[k]
      if k in (ALL_FEATURES + ['HGNC.symbol']): d[k] = entry[k]
    human_features.append(d)

  #remove all information aside from description and MCF7 label
  labels = {}
  for entry in human_labels:
    if entry['MCF7_BREAST'] != None:
      labels[entry['Description']] = 1 if float(entry['MCF7_BREAST']) <= -2 else 0

  #write to file
  human_data = {}
  sortedKeys = sorted(human_features[0].keys())
  sortedKeys.remove('HGNC.symbol')
  with open(out, 'w') as W:
    W.write(','.join(['name'] + sortedKeys + ['Essential'])+'\n')
    for entry in human_features:
      name = entry['HGNC.symbol']
      if name in labels:
        human_data[name] = {k: entry[k] for k in sortedKeys}
        human_data[name]['Essential'] = str(labels[name])
        s = ','.join([name] + [entry[k] for k in sortedKeys] + [str(labels[name])])
        W.write(s+'\n')

  return (human_data, sortedKeys)

if __name__ == '__main__':
  #(human_data, feature_set) = gen_human_data_file('../data/29Nov15_compiled_features_human.txt', '../data/essential_by_cancer_human.txt', '../data/human_data.csv')
  human_data = list(DictReader(open('../data/human_data.csv', 'r')))
  feature_set = [k for k in human_data[0].keys() if k != 'name' and k != 'Essential']

  #since there are only 5 features that the yeast and human data sets have in common, no need to filter the feature set
  #simply run the classifier on the human set and view the results

  args = classify.classify_args(features = feature_set,
                                classifiers=ALL_CLASSIFIERS,
                                write_to_log=True,
                                train_file="../data/testing_data_large.csv",
                                test_file="../data/human_data.csv",
                                classify=True)
  classify.main(args)

