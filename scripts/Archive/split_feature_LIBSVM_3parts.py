#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description='A script to split vector file into pfam/vfam feature, codon feature and 4-mer frequency feature')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
# parser.add_argument('-v', '--vector_number', help ='vector numbers for each feature to split, like 256,136 for 4-mer, 256,512 for 5-mer',required=True)

args = parser.parse_args()


# 30000: 30391

file_in_obj = open(args.input, 'r')
# vector_numbers = args.vector_number.split(",")

outfile_1 = open(args.input+".codon",'w')
outfile_2 = open(args.input+".4mers", 'w')
outfile_3 = open(args.input+".pvfam", 'w')

for line in file_in_obj:
    line = line.rstrip()
    fields = line.split()
    label = fields[0]




    file_1_line = label
    file_2_line = label
    file_3_line = label

    for vector in fields[1:]:
        feature_and_score = vector.split(":")
        if int(feature_and_score[0]) < 30000:
            file_3_line = file_3_line + ' ' + vector
        elif int(feature_and_score[0]) < 30256: # for codon, 30000 - 30255, 256 numbers
            file_1_line = file_1_line + ' ' + vector
        else:
            file_2_line = file_2_line + ' ' + vector

    outfile_1.write(file_1_line+'\n')
    outfile_2.write(file_2_line+'\n')
    if file_3_line != label:
        outfile_3.write(file_3_line+'\n')



