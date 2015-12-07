'''
Created on Dec 7, 2015

@author: Nicolas Metts

This module summarizes the number and percent of positive and negative labels
for a given data file
'''
import sys
import csv
from csv import DictReader

ESSENTIAL = 'Essential'


if __name__ == '__main__':
    data_file_name = sys.argv[1]
    data_file = open(data_file_name)
    data_reader = DictReader(data_file)
    pos_labels = 0
    neg_labels = 0
    for row in data_reader:
        if row[ESSENTIAL] == '0':
            neg_labels += 1
        else:
            pos_labels += 1
    data_file.close()
    data_summary_file = open('data/label_summary.csv', 'a')
    all_labels = pos_labels + neg_labels
    pos_percentage = float(pos_labels)/float(all_labels)
    neg_percentage = float(neg_labels)/float(all_labels)
    print "Summary for file: " + data_file_name
    print "\tTotal labels: " + str(all_labels)
    print "\tPositive labels: " + str(pos_labels)
    print "\tPositive percentage: " + str(pos_percentage)
    print "\tNegative percentage: " + str(neg_percentage)
    print "\tNegative labels: " + str(neg_labels)
    row = [data_file_name, pos_labels, pos_percentage, neg_labels, neg_percentage]
    data_summary_writer = csv.writer(data_summary_file)
    data_summary_writer.writerow(row)
    data_summary_file.close()