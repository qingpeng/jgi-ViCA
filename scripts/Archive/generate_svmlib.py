#!/usr/bin/env python
import argparse
import os
import math


parser = argparse.ArgumentParser(description='A script to convert vector(testing and training) file into LIBSVM format')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-o', '--output', help ='output file in LIBSVM format', required=True) 
parser.add_argument('-n', '--name', help ='name list of pfam/vfam', required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
file_out_obj = open(args.output, 'w')

name_file_obj = open(args.name, 'r')


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


name_list = []
for line in name_file_obj:
    line = line.rstrip()
    name_list.append(line.split()[1])
                
                
                
for line in file_in_obj:
    line =  line.rstrip()
    fields = line.split()
    newline = wordDic[fields[0]]
    
    if len(fields)>393:
        e_value = {}
        for field in fields[393:]:
            split_fields = field.split(":")
            label = name_list.index(split_fields[0])+1
            
            if not label in e_value:
                e_value[label] = float(split_fields[1])
            else:
                if float(split_fields[1])<e_value[label]:
                    e_value[label] = float(split_fields[1])
        sorted_key = sorted(e_value.keys())
        for label in sorted_key:
            if e_value[label] == 0.0:
                e_number = "1e-200"
            else:
                e_number = e_value[label]
#            print e_number
            field = str(label)+":"+str(math.log(float(e_number)))

            newline = newline + ' ' + field
            
    index = 30000
    for field in fields[1:393]:
        newline = newline + ' ' + str(index)+":"+field
        index = index + 1
        
    file_out_obj.write(newline+'\n')
    


