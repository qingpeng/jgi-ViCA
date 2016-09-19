#!/usr/bin/env python
import argparse
import os
import math

# for multiple label classification.... different taxonomical rank..

parser = argparse.ArgumentParser(description='A script to convert vector(testing and training) file into ARFF format')
parser.add_argument('-i', '--training', help ='Input vector file',required=True) 
parser.add_argument('-t', '--testing', help ='Input vector file',required=True) 

args = parser.parse_args()

file_training_obj = open(args.training, 'r')
file_testing_obj = open(args.testing, 'r')


file_training_out_obj = open(args.training+'.ARFF', 'w')
file_testing_out_obj = open(args.testing+'.ARFF', 'w')


# scan for 1st time, to get list of pfam id

key_set = set()

value_list_training = []
tax_list_training = []
for line in file_training_obj:
    line =  line.rstrip()
    fields = line.split('\t')
    newline = fields[0]
    value = {}
    if len(fields)>396: # 392 hmm+kmer vectors. + 4 columns in the beginning
        e_value = {}
        for field in fields[396:]: # for pfam and vfam
            split_fields = field.split(":")
            label = split_fields[0]
            
            if not label in e_value:
                e_value[label] = float(split_fields[1])
            else:
                if float(split_fields[1])<e_value[label]:
                    e_value[label] = float(split_fields[1])
        sorted_key = sorted(e_value.keys())
        for label in sorted_key: # covert to log .... for pfam/vfam e-value 
            if e_value[label] == 0.0:
                e_number = "1e-200"
            else:
                e_number = e_value[label]
            value[label] = math.log(float(e_number))
            key_set.add(label)
            newline = newline + ' ' + field
            
    index = 1
    for field in fields[4:396]:
        key_set.add(str(index))
        value[str(index)] = field
        index = index + 1
    value_list_training.append(value)
    taxonomy = fields[3].split()[2].split(',')[1].split('=')[1].split('/')
    tax_dict = {}
    for tax in taxonomy:
        if ':' in tax:
            fields = tax.split(':')
            tax_dict[fields[1]] = fields[0]
    tax = tax_dict['11']+'/'+tax_dict['16']+'/'+tax_dict['20']
    tax_list_training.append(tax)
    
    
value_list_testing = []
tax_list_testing = []

for line in file_testing_obj:
    line =  line.rstrip()
    fields = line.split('\t')
    newline = fields[0]
    value = {}
    if len(fields)>396: # 392 hmm+kmer vectors. + 4 columns in the beginning
        e_value = {}
        for field in fields[396:]: # for pfam and vfam
            split_fields = field.split(":")
            label = split_fields[0]
            
            if not label in e_value:
                e_value[label] = float(split_fields[1])
            else:
                if float(split_fields[1])<e_value[label]:
                    e_value[label] = float(split_fields[1])
        sorted_key = sorted(e_value.keys())
        for label in sorted_key: # covert to log .... for pfam/vfam e-value 
            if e_value[label] == 0.0:
                e_number = "1e-200"
            else:
                e_number = e_value[label]
            value[label] = math.log(float(e_number))
            key_set.add(label)
            newline = newline + ' ' + field
            
    index = 1
    for field in fields[4:396]:
        key_set.add(str(index))
        value[str(index)] = field
        index = index + 1
    value_list_testing.append(value)
    taxonomy = fields[3].split()[2].split(',')[1].split('=')[1].split('/')
    tax_dict = {}
    for tax in taxonomy:
        if ':' in tax:
            fields = tax.split(':')
            tax_dict[fields[1]] = fields[0]
    tax = tax_dict['11']+'/'+tax_dict['16']+'/'+tax_dict['20']
    tax_list_testing.append(tax)
    
    
    
line = "@RELATION Genelearn"
file_training_out_obj.write(line+'\n\n')
for key in key_set:
    line = "@ATTRIBUTE "+key+" numeric"
    file_training_out_obj.write(line+'\n')

tax_line = ','.join(list(set(tax_list_training)))
file_training_out_obj.write('@ATTRIBUTE class hierarchical '+tax_line+'\n')

file_training_out_obj.write('\n@DATA\n')

i = 0
for value in value_list_training:
    line = ''
    for key in key_set:
        if key in value:
            line = line + str(value[key])+','
        else:
            line = line +'?,'
    line = line  + tax_list_training[i]
    i += 1
    file_training_out_obj.write(line+'\n')



line = "@RELATION Genelearn"
file_testing_out_obj.write(line+'\n\n')
for key in key_set:
    line = "@ATTRIBUTE "+key+" numeric"
    file_testing_out_obj.write(line+'\n')

tax_line = ','.join(list(set(tax_list_training)))
file_testing_out_obj.write('@ATTRIBUTE class hierarchical '+tax_line+'\n')

file_testing_out_obj.write('\n@DATA\n')

i = 0
for value in value_list_testing:
    if tax_list_testing[i] in tax_list_training:
    
        line = ''
        for key in key_set:
            if key in value:
                line = line + str(value[key])+','
            else:
                line = line +'?,'
        line = line  + tax_list_testing[i]
    
        file_testing_out_obj.write(line+'\n')
    i += 1

