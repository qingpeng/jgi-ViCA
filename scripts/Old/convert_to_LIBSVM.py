#!/usr/bin/env python
import argparse
import os



parser = argparse.ArgumentParser(description='A script to convert vector file into LIBSVM format')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-o', '--output', help ='output file in LIBSVM format', required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
file_out_obj = open(args.output, 'w')


for line in file_in_obj:
    line =  line.rstrip()
    fields = line.split()
    
    if fields[0] == 'virus':
        newline = '1'
    else:
        newline = '0'
    index = 1
    for field in fields[1:]:
        newline = newline + ' ' + str(index)+":"+field
        index = index + 1
        
    file_out_obj.write(newline+'\n')


