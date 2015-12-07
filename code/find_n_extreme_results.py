'''
Created on Dec 7, 2015

@author: Nicolas Metts
'''
import sys
import csv
import numpy as np

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

PRECISION = 'Precision'
RECALL = 'Recall'
ACCURACY = 'Accuracy'

def __find_n_extreme_results(num_results, results_filename):
    results_file = open(results_filename)
    results_reader = csv.DictReader(results_file)
    field_names = results_reader.fieldnames
    accuracy_scores = []
    precision_scores = []
    recall_scores = []
    all_rows = []
    for row in results_reader:
        accuracy_scores.append(float(row[ACCURACY]))
        precision_scores.append(float(row[PRECISION]))
        recall_scores.append(float(row[RECALL]))
        all_rows.append(row)
    results_file.close()
    accuracy_scores = np.array(accuracy_scores)
    precision_scores = np.array(precision_scores)
    recall_scores = np.array(recall_scores)
    mean_accuracy = np.mean(accuracy_scores)
    mean_precision = np.mean(precision_scores)
    mean_recall = np.mean(recall_scores)
    sorted_accuracy_scores = np.argsort(accuracy_scores)
    sorted_recall_scores = np.argsort(recall_scores)
    sorted_precision_scores = np.argsort(precision_scores)
    top_n_accuracy_scores = sorted_accuracy_scores[len(accuracy_scores) - 1:
                                                   len(accuracy_scores) - num_results:
                                                   -1]
    bottom_n_accuracy_scores = sorted_accuracy_scores[0:num_results]
    top_n_precision_scores = sorted_precision_scores[len(precision_scores) - 1:
                                                   len(precision_scores) - num_results:
                                                   -1]
    bottom_n_precision_scores = sorted_precision_scores[0:num_results]
    top_n_recall_scores = sorted_recall_scores[len(recall_scores) - 1:
                                                   len(recall_scores) - num_results:
                                                   -1]
    bottom_n_recall_scores = sorted_recall_scores[0:num_results]
    top_results_file = open('logs/top_' + str(num_results) + '_results.csv', 'w')
    bottom_results_file = open('logs/bottom_' + str(num_results) +'_results.csv', 'w')
    top_results_writer = csv.DictWriter(top_results_file, delimiter=',', 
                                        fieldnames=field_names)
    bottom_results_writer = csv.DictWriter(bottom_results_file, delimiter=',',
                                           fieldnames=field_names)
    top_results_writer.writeheader()
    bottom_results_writer.writeheader()
    #top_results_writer.writerow(field_names)
    #bottom_results_writer.writerow(field_names)
    
    # Combine any overlapping top results
    all_top_results = set(np.concatenate((top_n_accuracy_scores,
                                          top_n_precision_scores,
                                          top_n_recall_scores)))
    all_bottom_results = set(np.concatenate((bottom_n_accuracy_scores,
                                          bottom_n_precision_scores,
                                          bottom_n_recall_scores)))
    for i in all_top_results:
        # Don't include top results if accuracy, precision, or recall is 0
        if (float(all_rows[i][ACCURACY]) > 0.0 and 
            float(all_rows[i][PRECISION]) > 0.0 and
            float(all_rows[i][RECALL]) > 0.0):
            top_results_writer.writerow(all_rows[i])
    for j in all_bottom_results:
        bottom_results_writer.writerow(all_rows[j])
    top_results_file.close()
    bottom_results_file.close()

if __name__ == '__main__':
    __find_n_extreme_results(int(sys.argv[1]), sys.argv[2])