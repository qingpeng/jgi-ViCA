#!/usr/bin/env python

import sys

file_in_obj = open(sys.argv[1], 'r')
file_virus_obj = open(sys.argv[2], 'w')
file_nonvirus_obj = open(sys.argv[3], 'w')

for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    label = fields[0][-1]
    if label == '0':
        file_nonvirus_obj.write(line+'\n')
    else:
        file_virus_obj.write(line+'\n')
