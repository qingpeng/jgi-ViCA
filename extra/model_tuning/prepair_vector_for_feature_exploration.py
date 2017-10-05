#!/usr/bin/env python
# split libsvm file into two files with label 1 and 0 , with feature score only
# to explore the pattern

import sys

file_in = sys.argv[1]
file_out_1 = file_in+'.label_1'
file_out_0 = file_in+'.label_0'

file_in_obj = open(file_in, 'r')

virus_out = open(file_out_1, 'w')
nonvirus_out = open(file_out_0, 'w')

for line in file_in_obj:
    line = line.rstrip()
    if line != '':
        fields = line.split()
        # print fields
        output = ''
        for field in fields[1:]:
            s = field.split(":")
            output = output+' '+s[1]
        if fields[0] == '0':
            nonvirus_out.write(output+'\n')
        else:
            virus_out.write(output+'\n')
