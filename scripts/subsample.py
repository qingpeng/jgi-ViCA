#!/usr/bin/env python


import sys

file_in_obj = open(sys.argv[1],'r')
file_out_obj = open(sys.argv[2],'w')


k = 0
for line in file_in_obj:
    if k == 10000: # or 20000
        k = 1
    else:
        k += 1
    if k <=200:
        file_out_obj.write(line)
