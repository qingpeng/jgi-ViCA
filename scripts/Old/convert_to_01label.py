#!/usr/bin/env python
import argparse
import os

# convert to 0/1 label

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
        label = '1'
    else:
        label = '0'
        
    newline = ' '.join([label]+fields[1:])
        
    file_out_obj.write(newline+'\n')


