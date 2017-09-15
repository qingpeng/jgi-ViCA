#!/usr/bin/env python


import sys
import random
file_in_obj = open(sys.argv[1], 'r')
virus_to_pick = int(sys.argv[2])
nonvirus_to_pick = int(sys.argv[3])

count_virus = 0
count_nonvirus = 0
for line in file_in_obj:
    if line[0] == '1':
        count_virus += 1
    else:
        count_nonvirus += 1

file_in_obj.close()

virus_ratio = virus_to_pick*1.0/count_virus
nonvirus_ratio = nonvirus_to_pick*1.0/count_nonvirus

random.seed(12345)

file_in_obj = open(sys.argv[1], 'r')
file_out_obj = open(sys.argv[4], 'w')


for line in file_in_obj:
    if line[0] == '1':
        random_num = random.random()
        if random_num < virus_ratio:
            file_out_obj.write(line)
    else:
        random_num = random.random()
        if random_num < nonvirus_ratio:
            file_out_obj.write(line)
