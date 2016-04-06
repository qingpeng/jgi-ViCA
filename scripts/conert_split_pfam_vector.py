#!/usr/bin/env python
import argparse
import os



parser = argparse.ArgumentParser(description='A script to convert vector(testing and training) file into LIBSVM format')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-o', '--output', help ='output file in LIBSVM format', required=True) 
parser.add_argument('-p', '--pfam', help ='output of pfam vector file in LIBSVM format', required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
file_out_obj = open(args.output, 'w')
file_pfam_obj = open(args.pfam, 'w')

wordDic = {
"Archaea" : "0",
"Bacteria" : "0",
"ssRNAPhage" : "1",
"ssRNAVirus" : "1",
"dsRNAPhage" : "1",
"dsRNAVirus" : "1",
"dsDNAPhage" : "1",
"dsDNAVirus" : "1",
"ssDNAPhage" : "1",
"ssDNAVirus" : "1",
"Retroviruses" : "1",
"Eukaryota" : "0",
"Mitochondrion" : "0",
"Chloroplast" : "0"}


for line in file_in_obj:
    line =  line.rstrip()
    fields = line.split()
    newline = wordDic[fields[0]]
    pfam_newline = newline
    index = 30000
    for field in fields[1:393]:
        newline = newline + ' ' + str(index)+":"+field
        index = index + 1
    
    if len(fields)>393:
        for field in fields[393:]:
            newline = newline + ' ' + field
            pfam_newline = pfam_newline + ' ' + field
        file_pfam_obj.write(pfam_newline+'\n')
    file_out_obj.write(newline+'\n')
    

