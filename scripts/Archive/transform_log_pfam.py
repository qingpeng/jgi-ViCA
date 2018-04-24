#!/usr/bin/env python
# logirize pfam e-value number
#


import sys
import math

file_in_obj = open(sys.argv[1], 'r')
file_out_obj = open(sys.argv[2], 'w')


for line in file_in_obj:
 #   print line
    line = line.rstrip()
    fields = line.split()
    newline = fields[0]
    for field in fields[1:]:
        f = field.split(":")
        # if e-value == 0.0, set it as 1e-200
        if int(f[0]) < 30000:
            if f[1] == "0.0":
                f1 = "1e-200"
            else:
                f1 = f[1]
            new_number = f[0] + ":" + str(math.log(float(f1)))
        else:
            new_number = field
        newline = newline + ' ' + new_number
        
    file_out_obj.write(newline+'\n')
        