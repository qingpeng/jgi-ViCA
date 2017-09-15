import sys
import random


file_vect = open(sys.argv[1], 'r')

id_set = set()

count = 0
name_list = []
length_list = []
for line in file_vect:
    line = line.rstrip()
    fields = line.split()
    name_list.append(fields[0])
    length_list.append(fields[1])
    if int(fields[1]) >= 2000:
        id_set.add(count)
    count += 1

file_prediction = open(sys.argv[2], 'r')

file_out = open(sys.argv[3], 'w')
count2 = 0
for line in file_prediction:
    line = line.rstrip()
    if count2 in id_set:
        to_print = name_list[count2]+' '+str(length_list[count2])+' '+line+'\n'
        file_out.write(to_print)
    count2 += 1
