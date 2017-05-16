#!/usr/bin/env python


import sys

file_in_obj = open(sys.argv[1],'r')
file_out1_obj = open(sys.argv[1]+'.virus','w')
file_out2_obj = open(sys.argv[1]+'.nonvirus','w')



for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    pairs = {}
    for field in fields[1:]:
        f = field.split(':')
        pairs[f[0]] = f[1]
    if '393' in pairs:
        to_write = pairs['393']+' '+pairs['394']+' '+pairs['395']+'\n'
        if fields[0] == '0':
            file_out2_obj.write(to_write)
        else:
            file_out1_obj.write(to_write)

