#!/usr/bin/env python
import argparse

#file_rank_obj = open("/global/projectb/scratch/qpzhang/Full_Training/Pfam/Vfam_run/all.vect.rank", 'r')
#file_out_obj = open("predictionAndLabels_family.txt.in_order", 'w')
#file_prediction_obj = open("predictionAndLabels_family.txt", 'r')
#file_family_index_obj = open("all.vect.family.index", 'r')

parser = argparse.ArgumentParser(description='A script to convert prediction to different level' )
parser.add_argument('-i', '--input', help ="prediction file to convert", required= True, default='predictionAndLabels_family.txt')
parser.add_argument('-f', '--family', help ="family.index file", required= True,default='all.vect.family.index')
parser.add_argument('-r', '--rank', help ="rank file with taxonomy level info", required= True,default='all.vect.rank')
parser.add_argument('-o', '--output', help ="output prediction file",required= True) 
args = parser.parse_args()


file_rank_obj = open(args.rank, 'r')
file_out_obj = open(args.output, 'w')
file_prediction_obj = open(args.input, 'r')
file_family_index_obj = open(args.family, 'r')


family_id = {}# 0,1,2 -> real tax-id
for line in file_family_index_obj:
    line = line.rstrip()
    fields = line.split()
    family_id[fields[0]] = fields[1]

superkindom_dict = {}


for line in file_rank_obj: # real tax-id family -> real order_id
    line = line.rstrip()
    tax_dict = {}
    taxonomy_list = line.split()[1].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
    if "0" in tax_dict and "16" in tax_dict:
            superkindom_dict[tax_dict['16']] = tax_dict['0']



for line in file_prediction_obj:# all family id from 0...]: 

    line = line.rstrip()
    fields = line[1:-2].split(",")
    newline = "(" + superkindom_dict[family_id[str(int(float(fields[0])))]] + \
              "," + superkindom_dict[family_id[str(int(float(fields[1])))]] + ")"
    file_out_obj.write(newline+'\n')


file_out_obj.close()


