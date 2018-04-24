#!/usr/bin/env python
import argparse
import os
parser = argparse.ArgumentParser(description='A script to deal with pfam vectors')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-o', '--output', help ='output vector file',required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
file_out_obj = open(args.output, 'w')

for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    e_value = {}
    for field in fields[1:]:
        s = field.split(":")
        if not int(s[0]) in e_value:
            e_value[int(s[0])] = float(s[1])
        else:
            if float(s[1])<e_value[int(s[0])]:
                e_value[int(s[0])]= float(s[1])
    sorted_key = sorted(e_value.keys())
    newline = fields[0]
    for key in sorted_key:
        newline = newline + ' ' + str(key)+":"+str(e_value[key])
    
    file_out_obj.write(newline+'\n')
    
        