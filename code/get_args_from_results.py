'''
Created on Dec 13, 2015

@author: Nicolas Metts
'''
import csv
from csv import DictReader
import sys

LARGE_ALL_FEATURES = ["Transcript.length", "Strand", "GC.content", "Enzyme",
                "SEG.low.complexity", "Transmembrane.domain", "Signal.peptide",
                "Coiled.coil", "Nucleus", "Mitochondria", "ER", "Cytoplasm",
                "Ribosome", "Interaction.partners", "Mw", "PI", "Protein.length",
                "Gravy.score", "Aromaticity.Score", "CAI", "Codon.Bias",
                "FOP.Score", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His",
                "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser",
                "Thr", "Val", "Trp", "Tyr", "Carbons", "Hydrogens", "Nitrogens",
                "Oxygens", "Sulphurs", "Instability.index", "Alphatic.index",
                "EC.number"]

SMALL_ALL_FEATURES = ["Mitochondria", "Cytoplasm", "ER", "Nucleus", "Vacuole",
                    "Other.location", "CAI", "NC", "GC.content",
                    "Protein.length", "Gravy.score", "DovEXPR", 
                    "BLAST.hits_in_yeast", "Interaction.partners", "Chromosome",
                    "Chr.position", "Intron", "Close.stop.ratio", 
                    "Rare.aa_ratio", "TM.helix", "in_how_many_of_5_proks", 
                    "in_how_many_of_6_close_yeast"]

def __get_args_from_results(results_file_name, args_file_name):
    results_file = open(results_file_name)
    results_reader = DictReader(results_file)
    args_file = open(args_file_name, 'w')
    for row in results_reader:
        args = []
        args.append("--classify")
        args.append("--data_file")
        args.append(row["Data_File"])
        args.append("--train_file")
        args.append(row["Training_Data"])
        args.append("--test_file")
        args.append(row["Testing_Data"])
        args.append("--classifiers")
        args.append(row["Classifier"])
        if row["Classifier"] == "svm":
            args.append("--kernel")
            args.append(row["Kernel"])
        if row["Scale"] == "TRUE":
            args.append("--scale")
        args.append("--features")
        num_features = 0
        if "large_yeast_data" in row["Data_File"]:
            for feature in LARGE_ALL_FEATURES:
                if row[feature] == "TRUE":
                    args.append(feature)
                    num_features += 1
        elif "small_yeast_data" in row["Data_File"]:
            for feature in SMALL_ALL_FEATURES:
                if row[feature] == "TRUE":
                    args.append(feature)
                    num_features += 1
        args_str = " ".join(args)
        args_file.write(args_str + "\n")
        print "Number of features: " + str(num_features)
    args_file.close()
    results_file.close()

if __name__ == '__main__':
    __get_args_from_results(sys.argv[1], sys.argv[2])