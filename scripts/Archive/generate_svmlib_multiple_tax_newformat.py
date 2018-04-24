#!/usr/bin/env python
import argparse
import os
import math

# for multiple label classification.... different taxonomical rank..


def fun_log(value):
    if value == 0.0:
        value = "1e-200"
    return math.log(float(value))
    




parser = argparse.ArgumentParser(description='A script to convert vector(testing and training) file into LIBSVM format')
parser.add_argument('-i', '--input', help ='Input vector file',required=True) 
parser.add_argument('-o', '--output', help ='output file in LIBSVM format', required=True) 
parser.add_argument('-n', '--feature', help ='file with list of feature and the id as the label in svmlib file', required=True) 

args = parser.parse_args()

file_in_obj = open(args.input, 'r')
file_out_obj = open(args.output, 'w')

feature_list_obj = open(args.feature, 'w')



                
list_label = []
            
for line in file_in_obj:
    line =  line.rstrip()
    fields = line.split('\t')
    
    tax_label = fields[0]
    vectors = fields[5]
    
    vectors_list = vectors.split(" ")
    
    vector_strings = []
    
    for vector in vectors_list:
        label = vector.split(":")[0]
        pre = vector.split("_")[0]
        value = float(vector.split(":")[1])
        if pre == "2" or pre == "3": # if pfam or vfam ,get log 
            value = fun_log(value)
        
        if not label in list_label:
            list_label.append(label)
        label_id = list_label.index(label)
        vector_strings.append(str(label_id)+":"+str(value))
    
    print_line = tax_label+' '+' '.join(vector_strings)
    
    file_out_obj.write(print_line+'\n')
    
        
id = 0
for label in list_label:
    print_line = str(id)+' '+label
    feature_list_obj.write(print_line+'\n')
    id += 1
    


