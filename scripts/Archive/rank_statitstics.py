#!/usr/bin/env python
import sys

file_in_obj = open(sys.argv[1], 'r')
file_out_obj = open(sys.argv[2], 'w')

count_full = 0
count_11_16 = 0
count_16_20 = 0
count_20_24 = 0

for line in file_in_obj:
    line = line.rstrip()
    tax_dict = {}
    taxonomy_list = line.split()[1].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
    if "11" in tax_dict and "16" in tax_dict and "20" in tax_dict and "24" in tax_dict:
        count_full += 1
    if "11" in tax_dict and "16" in tax_dict:
        count_11_16 += 1
    if "16" in tax_dict and "20" in tax_dict:
        count_16_20 += 1
    if "20" in tax_dict and "24" in tax_dict:
        count_20_24 += 1
    
file_out_obj.write("all_full:"+str(count_full)+'\n')
file_out_obj.write("11-16:"+str(count_11_16)+'\n')
file_out_obj.write("16-20:"+str(count_16_20)+'\n')
file_out_obj.write("20-24:"+str(count_20_24)+'\n')


