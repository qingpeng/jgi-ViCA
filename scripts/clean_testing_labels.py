#!/usr/bin/env python
# make sure all feature labels in testing data set is in training data set or model
#


import sys
file_training = open(sys.argv[1],'r')
file_testing = open(sys.argv[2],'r')

file_training_obj = open(sys.argv[1]+'.cleanup','w')
file_testing_obj = open(sys.argv[2] + '.cleanup','w')

for line in file_training:
    line = line.rstrip()
    fields = line.split()
    set_label_training = set()
    for vector in fields[1:]:
        label = vector.split(":")
        set_label_training.add(label[0])
        
for line in file_testing:
    line = line.rstrip()
    fields = line.split()
    set_label_test = set()
    for vector in fields[1:]:
        label = vector.split(":")
        set_label_test.add(label[0])
        

file_training = open(sys.argv[1],'r')
file_testing = open(sys.argv[2],'r')

for line in file_training:
    line = line.rstrip()
    fields = line.split()
    out = fields[0]
    for vector in fields[1:]:
        label = vector.split(":")
        if label[0] in set_label_training and label[0] in set_label_test:
            out = out + ' ' +vector
    file_training_obj.write(out+'\n')
    
    
    
for line in file_testing:
    line = line.rstrip()
    fields = line.split()
    out = fields[0]
    for vector in fields[1:]:
        label = vector.split(":")
        if label[0] in set_label_training and label[0] in set_label_test:
            out = out + ' ' +vector
    file_testing_obj.write(out+'\n')
    
    
