#!/usr/bin/env python

import sys


file_in_obj = open(sys.argv[1], 'r')
file_out_obj = open(sys.argv[2], 'w')

for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    newline = ' '.join(fields[1:])+'\n'
    file_out_obj.write(newline)