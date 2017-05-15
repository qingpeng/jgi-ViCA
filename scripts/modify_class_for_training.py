#!/usr/bin/env python

import sys
import random

file_old_vector_obj = open(sys.argv[1], 'r')
file_new_vector_obj = open(sys.argv[2], 'r')
file_out_obj = open(sys.argv[3], 'w')

class_list = []

for line in file_old_vector_obj:
    fields = line.split()
    class_list.append(fields[0])

i = 0
for line in file_new_vector_obj:
    fields = line.split()
    newline = class_list[i]+' '+' '.join(fields[1:])
    file_out_obj.write(newline+'\n')

    i += 1


