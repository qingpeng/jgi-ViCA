#!/usr/bin/env python

import sys


num = int(sys.argv[2])
cntdict ={}
with open(sys.argv[1], 'r') as file:
    for line in file:
        ln = line.split('\t')
        type = ln[0]
        if type in cntdict:
            if cntdict[type] <= num:
                cntdict[type] = cntdict[type] + 1 
                print(line)
            else:
                continue
        else:        
            cntdict[type] = 1