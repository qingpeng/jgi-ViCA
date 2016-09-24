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
        else: # if there are <2 sub nodes, return nothing
            dictionary_rank_testing[key] = []
    return dictionary_rank_testing

parser = argparse.ArgumentParser(description='A script to split vector file into training and testing by different tax rank')
parser.add_argument('-v', '--vector', help ='full vector file',required=True)
#parser.add_argument('-r', '--rank', help ='rank file',required=True) 
#parser.add_argument('-o', '--outfile', help ='output file', required=True) 

args = parser.parse_args()


file_vector_obj = open(args.vector, 'r')
#file_rank_obj = open(args.rank, 'r')

file_order_training_obj = open(args.vector+'.order.training', 'w')
file_family_training_obj = open(args.vector+'.family.training', 'w')
file_genus_training_obj = open(args.vector+'.genus.training', 'w')

file_order_testing_obj = open(args.vector+'.order.testing', 'w')
file_family_testing_obj = open(args.vector+'.family.testing', 'w')
file_genus_testing_obj = open(args.vector+'.genus.testing', 'w')

file_order_index_obj = open(args.vector + '.order.index', 'w')
file_family_index_obj = open(args.vector + '.family.index', 'w')
file_genus_index_obj = open(args.vector + '.genus.index', 'w')

order = {}

family = {}

genus = {}



# 11 - order
# 16 - family
# 20 - genus
# 24 - species
#The
#          "taxonomy" option outputs the complete taxonomic path in the
#          format: "taxonId0:rankId0/..." (note that not all nodes have
#          rankIds defined).
         


for line in file_vector_obj:
    line = line.rstrip()
    tax_dict = {}
#    print line
    taxonomy_list = line.split('\t')[2].split(',')[1].split("=")[1].split("/")
#    print taxonomy_list
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
            
print "1st time scanning vector file done!\n"
file_vector_obj.close()

file_vector_obj = open(args.vector, 'r')


# get the list of family id for testing for each order id
order_testing = sample_testing(order)

# the same for each family id
family_testing = sample_testing(family)
genus_testing = sample_testing(genus)

#exit()
# scan original vector file, and put the lines into different files for testing and training for different family, order, gener

# the list of order, family, genus ids... in output files, the labels are integer... from 0,,,k-1. 
# so need to convert to integer from real tax id, the tax id and integer relation can be
# retrieved from the index output file
#
order_list = []
family_list = []
genus_list = []



for line in file_vector_obj:
    line = line.rstrip()
    fields = line.split('\t')
    tax_dict = {}
#    print fields
    taxonomy_list = fields[2].split(',')[1].split("=")[1].split("/")
    for taxonomy in taxonomy_list:
        groups = taxonomy.split(":")
        if groups[1] != "":
            tax_dict[groups[1]] = groups[0]
 #   print tax_dict
    if "11" in tax_dict and "16" in tax_dict and "20" in tax_dict and "24" in tax_dict:
    # must have all 4 tax id in 4 levels
    
    
    
        if tax_dict["16"] in order_testing[tax_dict["11"]]: 
        # if the family of this  record is picked as testing for the order it belongs to
            try:
                order_id_num = order_list.index(tax_dict["11"])
            except ValueError:
                order_id_num = len(order_list)
                order_list.append(tax_dict["11"])
            file_order_testing_obj.write(str(order_id_num)+'\t'+tax_dict["16"]+'\t'+line+'\n')
        else: # if this record is for training
            if order_testing[tax_dict["11"]] != []: # only if there are >=2 sub nodes

                try:
                    order_id_num = order_list.index(tax_dict["11"])
                except ValueError:
                    order_id_num = len(order_list)
                    order_list.append(tax_dict["11"])
                
                
                file_order_training_obj.write(str(order_id_num)+'\t'+tax_dict["16"]+'\t'+line+'\n')
        # output for order level...
        # output format : order_id_number family_id ...
        
                
        if tax_dict["20"] in family_testing[tax_dict["16"]]:
        
            try:
                family_id_num = family_list.index(tax_dict["16"])
            except ValueError:
                family_id_num = len(family_list)
                family_list.append(tax_dict["16"])
                
            file_family_testing_obj.write(str(family_id_num)+'\t'+tax_dict["20"]+'\t'+line+'\n')
        else:
            if family_testing[tax_dict["16"]] != []:

                try:
                    family_id_num = family_list.index(tax_dict["16"])
                except ValueError:
                    family_id_num = len(family_list)
                    family_list.append(tax_dict["16"])
                

                file_family_training_obj.write(str(family_id_num)+'\t'+tax_dict["20"]+'\t'+line+'\n')
                
                
        if tax_dict["24"] in genus_testing[tax_dict["20"]]:

            try:
                genus_id_num = genus_list.index(tax_dict["20"])
            except ValueError:
                genus_id_num = len(genus_list)
                genus_list.append(tax_dict["20"])
                
                
            file_genus_testing_obj.write(str(genus_id_num)+'\t'+tax_dict["24"]+'\t'+line+'\n')
        else:
            if genus_testing[tax_dict["20"]] != []:

                try:
                    genus_id_num = genus_list.index(tax_dict["20"])
                except ValueError:
                    genus_id_num = len(genus_list)
                    genus_list.append(tax_dict["20"])
                
                
                file_genus_training_obj.write(str(genus_id_num)+'\t'+tax_dict["24"]+'\t'+line+'\n')


for i in range(len(order_list)):
    file_order_index_obj.write(str(i) + ' ' + order_list[i] + '\n')

for i in range(len(family_list)):
    file_family_index_obj.write(str(i) + ' ' + family_list[i] + '\n')

for i in range(len(genus_list)):
    file_genus_index_obj.write(str(i) + ' ' + genus_list[i] + '\n')




