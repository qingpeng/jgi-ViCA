#!/usr/bin/env python
import argparse
import random



# split vector into testing data and training data:
# criteria:
# 1. all four ranks exist
# 2. at least 2 subgroups.
# 3. if <=5 subgroups, 1 as testing, other as training
# 4. if >5 subgroups, ,20% as testing +1 , other as training
# 5.

def sample_testing(dictionary_rank):
# return dictioanryof list with the testing taxid for each tax id in that rank
    dictionary_rank_testing = {}
    for key in dictionary_rank:
        length = len(dictionary_rank[key])
        if length >= 2: # only if there are >=2 sub nodes in any taxon id
            if length <=5:
                number_to_sample = 1
            else:
                number_to_sample = length/5 + 1
            list_of_subnodes = list(dictionary_rank[key]) # convert from set to list to shuffle
            random.shuffle(list_of_subnodes)
            dictionary_rank_testing[key] =  list_of_subnodes[:number_to_sample]
        else:
            dictionary_rank_testing[key] = []
    return dictionary_rank_testing

parser = argparse.ArgumentParser(description='A script to split vector file into training and testing by different tax rank')
parser.add_argument('-v', '--vector', help ='full vector file',required=True)
parser.add_argument('-r', '--rank', help ='rank file',required=True) 
#parser.add_argument('-o', '--outfile', help ='output file', required=True) 

args = parser.parse_args()


file_vector_obj = open(args.vector, 'r')
file_rank_obj = open(args.rank, 'r')

file_order_training_obj = open(args.vector+'.order_training', 'w')
file_family_training_obj = open(args.vector+'.family_training', 'w')
file_genus_training_obj = open(args.vector+'.genus_training', 'w')

file_order_testing_obj = open(args.vector+'.order_testing', 'w')
file_family_testing_obj = open(args.vector+'.family_testing', 'w')
file_genus_testing_obj = open(args.vector+'.genus_testing', 'w')


order = {}

family = {}

genus = {}



for line in file_rank_obj:
    line = line.rstrip()
    tax_dict = {}
    taxonomy_list = line.split()[1].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
    if "11" in tax_dict and "16" in tax_dict and "20" in tax_dict and "24" in tax_dict:

        try:
            order[tax_dict["11"]].add(tax_dict["16"])
        except KeyError:
            order[tax_dict["11"]] = set([tax_dict["16"]]) # use set to remove duplicate
        try:
            family[tax_dict["16"]].add(tax_dict["20"])
        except KeyError:
            family[tax_dict["16"]] = set([tax_dict["20"]])
        try:
            genus[tax_dict["20"]].add(tax_dict["24"])
        except KeyError:
            genus[tax_dict["20"]] = set([tax_dict["24"]])
            
print "loading rank file done!\n"


order_testing = sample_testing(order)
#print family
family_testing = sample_testing(family)
genus_testing = sample_testing(genus)

#exit()
# scan original vector file, and put the lines into different files for testing and training for different family, order, gener
for line in file_vector_obj:
    line = line.rstrip()
    fields = line.split()
    tax_dict = {}
#    print fields
    taxonomy_list = fields[3].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
 #   print tax_dict
    if "11" in tax_dict and "16" in tax_dict and "20" in tax_dict and "24" in tax_dict:
        if tax_dict["16"] in order_testing[tax_dict["11"]]:
            file_order_testing_obj.write(tax_dict["11"]+'\t'+tax_dict["16"]+'\t'+line+'\n')
        else:
            if order_testing[tax_dict["11"]] != []:
                file_order_training_obj.write(tax_dict["11"]+'\t'+tax_dict["16"]+'\t'+line+'\n')
        if tax_dict["20"] in family_testing[tax_dict["16"]]:
            file_family_testing_obj.write(tax_dict["16"]+'\t'+tax_dict["20"]+'\t'+line+'\n')
        else:
            if family_testing[tax_dict["16"]] != []:
                file_family_training_obj.write(tax_dict["16"]+'\t'+tax_dict["20"]+'\t'+line+'\n')
        if tax_dict["24"] in genus_testing[tax_dict["20"]]:
            file_genus_testing_obj.write(tax_dict["20"]+'\t'+tax_dict["24"]+'\t'+line+'\n')
        else:
            if genus_testing[tax_dict["20"]] != []:
                file_genus_training_obj.write(tax_dict["20"]+'\t'+tax_dict["24"]+'\t'+line+'\n')

