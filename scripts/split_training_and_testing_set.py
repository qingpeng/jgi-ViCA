import sys
import random


file_in_obj = open(sys.argv[1], 'r')

ratio = float(sys.argv[2])

file_training = open(sys.argv[3], 'w')
file_testing = open(sys.argv[4], 'w')

for line in file_in_obj:
    if random.random() <= ratio:
        file_training.write(line)
    else:
        file_testing.write(line)


