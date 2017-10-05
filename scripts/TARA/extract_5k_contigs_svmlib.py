import sys
import random


file_vect = open(sys.argv[1], 'r')  # vect file

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

file_libsvm = open(sys.argv[2], 'r')  # libsvm file

file_out = open(sys.argv[3], 'w')
count2 = 0
for line in file_libsvm:
    if count2 in id_set:
        file_out.write(line)
    count2 += 1
