#!/usr/bin/env python


import sys
import random

to_pick = 150600.0/5076200
random.seed(1)

file_in_obj = open(sys.argv[1], 'r')

file_out1_obj = open(sys.argv[1]+'.1x_nonvirus', 'w')
file_out2_obj = open(sys.argv[1]+'.2x_nonvirus', 'w')
file_out3_obj = open(sys.argv[1]+'.3x_nonvirus', 'w')

for line in file_in_obj:
    if line[0] == '1':
        file_out1_obj.write(line)
        file_out2_obj.write(line)
        file_out3_obj.write(line)

    else:
        random_num = random.random()
        if random_num < to_pick:
            file_out1_obj.write(line)
        if random_num < to_pick*2:
            file_out2_obj.write(line)
        if random_num < to_pick*3:
            file_out3_obj.write(line)
