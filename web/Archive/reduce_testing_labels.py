#!/usr/bin/env python
# make sure all feature labels in testing data set is in training data set or model
#


import sys
file_training = open(sys.argv[1],'r')
file_testing = open(sys.argv[2],'r')

file_testing_obj = open(sys.argv[2] + '.cleanup','w')

max_label = 0
for line in file_training:
    line = line.rstrip()
    fields = line.split()

    for vector in fields[1:]:
        label = vector.split(":")
        if int(label[0]) > max_label:
            max_label = int(label[0])

print max_label
        

    
for line in file_testing:
    line = line.rstrip()
    fields = line.split()
    out = fields[0]
    for vector in fields[1:]:
        label = vector.split(":")
        if int(label[0]) <= max_label:
            out = out + ' ' +vector
    file_testing_obj.write(out+'\n')
    
    
