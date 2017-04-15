#!/usr/bin/env python

import sys
import random

file_in_obj = open(sys.argv[1], 'r')
row_num = int(sys.argv[2])
pick_num = int(sys.argv[3])

row_pick_set = set(random.sample(xrange(row_num), pick_num))

file_out_obj = open(sys.argv[4], 'w')

n = 0

for line in file_in_obj:
    if n in row_pick_set:
        line = line.rstrip()
        fields = line.split()
        seg_id = fields[0][0:-1]
        label = fields[0][-1]
        newline = seg_id+' '+label+' '+' '.join(fields[1:])+'\n'
        file_out_obj.write(newline)

    n = n+1

