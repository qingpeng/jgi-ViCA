#!/usr/bin/env python

import sys
file_in = open(sys.argv[1],'r')

for line in file_in:
    line = line.rstrip()
    fields = line.split()
    out = fields[0]
    dict_vectors = {}
    for vector in fields[1:]:
        label = vector.split(":")
        dict_vectors[int(label[0])+1] = label[1]
    
    sorted_key = sorted(dict_vectors.keys())
    for key in sorted_key:
        out = out + ' '+ str(key)+':'+dict_vectors[key]
    
    print out
