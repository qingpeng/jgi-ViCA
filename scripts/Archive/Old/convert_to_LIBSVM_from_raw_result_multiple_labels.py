#!/usr/bin/env python
import argparse
import os



parser = argparse.ArgumentParser(description='A script to convert vector(testing and training) file into LIBSVM format for multiple classification')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-o', '--output', help ='output file in LIBSVM format', required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
file_out_obj = open(args.output, 'w')

wordDic = {
"Archaea" : "0",
"Bacteria" : "1",
"ssRNAPhage" : "2",
"ssRNAVirus" : "3",
"dsRNAPhage" : "4",
"dsRNAVirus" : "5",
"dsDNAPhage" : "6",
"dsDNAVirus" : "7",
"ssDNAPhage" : "8",
"ssDNAVirus" : "9",
"Retroviruses" : "10",
"Eukaryota" : "11",
"Mitochondrion" : "12",
"Chloroplast" : "13"}


for line in file_in_obj:
    line =  line.rstrip()
    fields = line.split()
    newline = wordDic[fields[0]]

    index = 1
    for field in fields[1:]:
        newline = newline + ' ' + str(index)+":"+field
        index = index + 1
        
    file_out_obj.write(newline+'\n')


