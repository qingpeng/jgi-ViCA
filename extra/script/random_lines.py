import sys
import random


file_in = sys.argv[1]
line_total = int(sys.argv[2])
line_to_pick = int(sys.argv[3])


to_pick_range = float(line_to_pick)/line_total


file_in_obj = open(file_in,'r')


for line in file_in_obj:
    line = line.rstrip()
    if random.random() <= to_pick_range:
        print line

