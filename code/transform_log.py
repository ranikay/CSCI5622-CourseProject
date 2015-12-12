'''
Created on Dec 12, 2015

@author: Nicolas Metts

A class to add a hash entry to a CSV log that doesn't have one

'''
import csv
from csv import DictReader
import sys

LEAVE_OUT_OF_HASH = ['Train_Size', 'Test_Size', 'Precision', 'Recall',
                     'Accuracy', 'True_Count', 'Actual_Count']

def __transform_log(file_name):
    old_file = open(file_name)
    old_file_reader = DictReader(old_file)
    old_header = old_file_reader.fieldnames
    if 'hash' in old_header:
        print file_name + " already has hash"
        return
    else:
        print "Adding hash to file: " + file_name
    old_contents = list(old_file_reader)
    old_file.close()
    with open(file_name, 'w') as new_file:
        new_writer = csv.writer(new_file)
        new_header = ['hash']
        new_header += old_header
        new_writer.writerow(new_header)
        for row in old_contents:
            row_contents = []
            for h in old_header:
                if h not in LEAVE_OUT_OF_HASH:
                    row_contents.append(row[h])
            row_str = "_".join([str(x) for x in row_contents])
            predict_hash = hash(row_str)
            new_row = [predict_hash]
            new_row += [row[h] for h in old_header]
            new_writer.writerow(new_row)
    
    print "Done modifying file: " + file_name
            
if __name__ == '__main__':
    file_name = sys.argv[1]
    __transform_log(file_name)