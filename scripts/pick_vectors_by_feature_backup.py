#!/usr/bin/env python


import argparse

parser = argparse.ArgumentParser(description='A script to pick vectors by feature id')
parser.add_argument('-d', '--index', help ='feature_index file',required=True)
parser.add_argument('-f', '--feature', help ='feature id to pick(for now, only 1 accepted, to remove)',required=True)
parser.add_argument('-i', '--input', help ='input vector file',required=True) 
parser.add_argument('-o', '--outfile', help ='output file', required=True) 


args = parser.parse_args()


file_index_obj = open(args.index, 'r')
feature = args.feature
file_in_obj = open(args.input, 'r')
file_output_obj = open(args.outfile, 'w')

index_set = set()
for line in file_index_obj:
    line =  line.rstrip()
    fields = line.split()
    if fields[1].split("_")[0] == feature:
        index_set.add(fields[0])
        
#print index_set
for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    newline = fields[0]
    for vector in fields[1:]:
        items = vector.split(":")
        if items[0] not in index_set:
            newline = newline + ' ' + vector
    file_output_obj.write(newline+'\n')
    


