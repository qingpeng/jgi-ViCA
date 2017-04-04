#!/usr/bin/env python

import sys

file_segment = sys.argv[1]
file_order_obj = open(sys.argv[2], 'r')
file_family_obj = open(sys.argv[3], 'r')
file_genus_obj = open(sys.argv[4], 'r')

out_order_obj = open(sys.argv[2]+'.valid.m8', 'w')
out_family_obj = open(sys.argv[3]+'.valid.m8', 'w')
out_genus_obj = open(sys.argv[4]+'.valid.m8', 'w')




file_segment_obj = open(file_segment, 'r')

family_set = set()
order_set = set()
genus_set = set()

count = 1
for line in file_segment_obj:
    if line[0] == ">":
        line = line.rstrip()
        print line
        taxonomy_list1 = line.split()
        print taxonomy_list1
        taxonomy_list = taxonomy_list1[2].split(',')[1].split("=")[1].split("/")
        #    print taxonomy_list
        tax_dict = {}
        for taxonomy in taxonomy_list:
            groups = taxonomy.split(":")
            if groups[1] != "":
                tax_dict[groups[1]] = groups[0]
                # 16 family 11 - order, 20 - genus

        if "11" in tax_dict:
            order_set.add(count)
        if "16" in tax_dict:
            family_set.add(count)
        if "20" in tax_dict:
            genus_set.add(count)

        count += 1

print len(order_set)
print len(family_set)
print len(genus_set)


for line in file_order_obj:
    line = line.rstrip()
    name = line.split()[0]
    if int(name) in order_set:
        out_order_obj.write(line+'\n')

for line in file_family_obj:
    line = line.rstrip()
    name = line.split()[0]
    if int(name) in family_set:
        out_family_obj.write(line+'\n')

for line in file_genus_obj:
    line = line.rstrip()
    name = line.split()[0]
    if int(name) in genus_set:
        out_genus_obj.write(line+'\n')

