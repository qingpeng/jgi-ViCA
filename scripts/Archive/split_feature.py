#!/usr/bin/env python
import argparse
import os

parser = argparse.ArgumentParser(description='A script to split vector file into codon feature and 4-mer frequency feature')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-v', '--vector_number', help ='vector numbers for each feature to split, like 256,136 256,512 for 5-mer',required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
vector_numbers = args.vector_number.split(",")

outfile_1 = open(args.input+".codon",'w')
outfile_2 = open(args.input+".4mers", 'w')


for line in file_in_obj:
        line =  line.rstrip()
        fields = line.split()
    
	file_1_line = fields[0]
	index = 1
	for i in range(1,1+int(vector_numbers[0])):
		file_1_line = file_1_line+' '+str(index)+":"+fields[i].split(":")[1]
		index = index+1
	outfile_1.write(file_1_line+'\n')
	
	
	file_2_line = fields[0]
	index = 1
	for i in range(1+int(vector_numbers[0]),1+int(vector_numbers[0])+int(vector_numbers[1])):
		file_2_line = file_2_line+' '+str(index)+":"+fields[i].split(":")[1]
		index = index+1
	outfile_2.write(file_2_line+'\n')
	

