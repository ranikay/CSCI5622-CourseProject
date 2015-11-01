import argparse
from numpy import random

#given a data file, split its rows into training and testing files
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Machine learning classifier options')
  parser.add_argument('--pTrain', type=float, default=0.15,
                      help="Separate data into training and testing sets randomly based on the probability of an example being in the training set")
  parser.add_argument('--trainFile', type=str, default='../data/training_data.csv',
                      help="Name of the resutling training file")
  parser.add_argument('--testFile', type=str, default='../data/testing_data.csv',
                      help="Name of the resulting testing file")
  parser.add_argument('--dataFile', type=str, default='../data/cerevisiae_compiled_features.csv',
                      help="Name of data file to split")
  args = parser.parse_args()

  with open(args.dataFile, 'rU') as data, open(args.trainFile, 'w') as train, open(args.testFile, 'w') as test:
    #assume the first line is column names
    cols = data.readline()
    train.write(cols)
    test.write(cols)
    for line in data:
      if random.rand() > args.pTrain: test.write(line)
      else: train.write(line)

