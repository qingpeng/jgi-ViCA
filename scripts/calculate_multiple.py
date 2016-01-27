#!/usr/bin/env python

import sys

file_in = open(sys.argv[1],'r')
file_out = open('multiple.out','w')

file_in.readline()
file_in.readline()
file_in.readline()
virus = ["ssRNAPhage","ssRNAVirus", "dsRNAPhage","dsRNAVirus","dsDNAPhage", "dsDNAVirus", "ssDNAPhage", "ssDNAVirus", "Retroviruses"]
fp=0
fn=0
tp=0
tn=0

for line in file_in:
    line = line.rstrip()
    fields = line.split()
    print fields
    if fields[1] in virus:
        if fields[3] in virus :
            tp = tp+int(fields[5])
        else:
            fn = fn+int(fields[5])
    else:
        if fields[3] in virus:
            fp = fp+int(fields[5])
        else:
            tn = tn+int(fields[5])

print "tp:"+str(tp)
print "fn:"+str(fn)
print "fp:"+str(fp)
print "tn:"+str(tn)

