# Integrating heterogeneous data sources to predict gene essentiality in Saccharomyces cerevisiae
Course project for CU Boulder CSCI 5622: Machine Learning with Dr. Jordan Boyd-Graber

**Team members**: Nicolas Metts, Matthew Pennington, Rani Schwindt and Carter Tillquist

### Introduction
An "essential" gene is one which, when absent/deleted, confers a lethal phenotype. We propose that gene essentiality can be predicted using a weighted combination of features. Here we use data sets with genes from *Saccharomyces cerevisiae*, a well characterized yeast species. 

### Seringhaus data set (data/cerevisiae_compiled_features.csv)

In 2006, Seringhaus _et al_ used 14 biological features to train a classifier for predicting essential genes in _S. cerevisiae_ and a related organism _S. miktae_. On 4,648 genes in _S. cerevisiae_, the classifier resulted in a precision TP/TP+FP = 0.69 and recall TP/TP+FN = 0.091$. The classifier used was an average of 7 different classifiers, including logistic regression, Naive Bayes, and AdaBoost. 

The data from this paper was provided by the lab [here](http://www.gersteinlab.org/proj/predess/data/Scerevisiae/Compiled/cerevisiae_ALL_noNaN.csv), however, it only includes a complete feature set for 3,500 genes. 

Label = SGD_ess (1 = essential, 0 = nonessential)

#### Seringhaus features for 3,500 genes

|       Feature        |      Description                                         |       Raw data format        |
|----------------------|----------------------------------------------------------|------------------------------|
| Mitochondria         | Does the protein localize to the mitochondria (predicted)|            Binary            | 
| Cytoplasm            | Does the protein localize to the cytoplasm (predicted)   |            Binary            | 
| ER                   | Does the protein localize to the ER (predicted)          |            Binary            | 
| Nucleus              | Does the protein localize to the nucleus (predicted)     |            Binary            | 
| Vacuole              | Does the protein localize to the vacuole (predicted)     |            Binary            | 
| Other                | Does the protein localize somewhere else (predicted)     |            Binary            | 
| CAI                  | [Codon adaptation index](http://en.wikipedia.org/wiki/Codon_Adaptation_Index)|  0 to 1  | 
| Nc                   | [Effective number of codons](http://en.wikipedia.org/wiki/Effective_number_of_codons)| Integer |
| GC                   | [GC content](http://www.pnas.org/content/111/39/E4096.long)|   0 to 1        |
| L_aa                 | Number of amino acids in protein (predicted)             |            Integer           |
| Gravy                | Hydrophobicity (positive) or hydrophilicity (negative)   |           -inf to inf        |
| DovEXPR              | Unknown                                                  |          0 to inf            |
| BLAST_hits_in_yeast  | Number of related genes in yeast (BLAST similarity)      |            Integer           |
| INTXN_partners       | Number of protein interaction partners                   |            Integer           |
| Chromosome           | Which chromosome (yeast have 16) is the gene on          |            Integer           |
| Chr_position         | Where does the gene start relative to the whole chromosome|            0 to 1           |
| Intron               | Unknown                                                  |            Binary            |
| CLOSE_STOP_RATIO     | % of codons one third-base away from stop codon          |            0 to 1            |
| RARE_AA_RATIO        | % of rare amino acids in translated ORF                  |            0 to 1            |
| TM_HELIX             | Number of transmembrane helices (predicted)              |            Integer           |
| In_how_many_of_5_proks_BLAST| Number of related genes in 5 prokaryotes (BLAST similarity)|            Integer  |
| In_how_many_of_6_close_yeast_BLAST| Number of related genes in 6 yeast species (BLAST similarity)|     Integer |

### Rani's data set in progress (data/07Nov15_all_genes_features.txt)

Compiling features for all 5,799 genes in _S. cerevisiae_. Label = Essential (1 = essential, 0 = nonessential). Note that the localization features in this data set are not predicted, but were experimentally determined.

#### Features for 5,799 genes

|       Feature        |      Description                                         |       Raw data format        |
|----------------------|----------------------------------------------------------|------------------------------|
| Transcript length    | Length of the transcribe gene including UTRs             |            Integer           |
| Strand               | Whether gene is on the positive DNA strand (1) or negative (-1)|            1 or -1     |
| GC                   | [GC content](http://www.pnas.org/content/111/39/E4096.long)|           0 to 1           |
| Enzyme               | Does the protein have enzymatic activity                 |             0 to 1           |
| SEG.low.complexity   | Predicted to have [low-complexity regions](http://mendel.imp.ac.at/METHODS/seg.server.html)|0 to 1|
| Transmembrane.domain | Does the protein have a [transmembrane domain](http://www.cbs.dtu.dk/services/TMHMM/)|0 to 1|
| Signal.peptide       | Does the protein have a [signal peptide](http://www.cbs.dtu.dk/services/SignalP/)|0 to 1|
| Coiled.coil          | Does the protein have a [coiled coil](http://www.ch.embnet.org/software/COILS_form.html)|0 to 1|
| Nucleus              | Does the protein localize to the nucleus                 |            Binary            |
| Mitochondria         | Does the protein localize to the mitochondria            |            Binary            | 
| ER                   | Does the protein localize to the ER                      |            Binary            | 
| Cytoplasm            | Does the protein localize to the cytoplasm               |            Binary            | 
| Ribosome             | Does the protein localize to the ribosome                |            Binary            | 

TO DO: add Gene Ontology terms about protein and gene function

