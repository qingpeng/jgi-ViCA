#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
import sys
import numpy


parser = argparse.ArgumentParser( description='A script to select taxonomic groups with particular attributes from a reftree feature vector database and split the data into test and train datasets' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-o', '--output', help ="An output vector file", type=argparse.FileType('w'), default='-')
parser.add_argument('-t', '--testfile', help ="An output vector file", type=argparse.FileType('w'), default='testdata.txt')
parser.add_argument('-r', '--reftree', help ="Reftree database directory location") 
args = parser.parse_args()

leveldict ={ "no_rank": None, 
"superkingdom" : 0,
"kingdom" : 1,
"subkingdom" : 2,
"superphylum" : 3,
"phylum" : 4,
"subphylum" : 5,
"superclass" : 6,
"class" : 7,
"subclass" : 8,
"infraclass" : 9,
"superorder" : 10,
"order" : 11,
"suborder" : 12,
"infraorder" : 13,
"parvorder" : 14,
"superfamily" : 15,
"family" : 16,
"subfamily" : 17,
"tribe" : 18,
"subtribe" : 19,
"genus" : 20,
"subgenus" : 21,
"species_group" : 22,
"species_subgroup" : 23,
"species" : 24,
"subspecies" : 25,
"varietas" : 26,
"forma" : 27
}

# Read the configuration file
config = json.load(open(args.config, 'r'))["training_data_formatter"]

# Variables
dbdir,dbname = os.path.split(os.path.abspath(args.reftree))
testlist = []
trainlist =[]

# Functions

#todo add level and attribute filtering once Ed adds it to Reftree
def get_reftree_data(reftreeopts):
    """A function to run RefTree"""
    p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE)
    reftreeout, reftreeerr= p0.communicate()
    assert p0.returncode == 0, "RefTree returned an error while searching the tree with taxid"
    return reftreeout

def organelle(line):
    """Returns the organelle value of a line"""
    attributedict = {}
    lv = line.strip().split('\t')
    desc = lv[2].split(' ')[1].split(',')
    for item in desc:
        name, var = line.partition("=")[::2]
        attributedict[name.strip()] = var
        if "organelle" in attributedict:
            return attributedict["organelle"]
        else:
            return None


def split_data(key, line, trainhandle, testhandle, trainlist, testlist, testprop):
    """takes a vector and assigns it to the test or train datasets updating the taxid lists """
    lv = line.strip().split('\t')
    taxid = lv[0]
    simple = [key] + lv[3:]
    if taxid in trainlist:
        trainhandle.write('\t'.join(simple) + '\n')
        return (trainlist, testlist)
    if taxid in testlist:
        testhandle.write('\t'.join(simple) + '\n')
        return (trainlist, testlist)
    else:
        selection = numpy.random.choice(["trainlist","testlist"], p=[1-testprop,testprop])
        if selection == "trainlist":
            trainlist.append(taxid)
            trainhandle.write(simple)
            return (trainlist, testlist)
        else:
            testlist.append(taxid)
            testhandle.write(simple)
            return (trainlist, testlist)

# MAIN

# If a taxonomic rank is specified return the data for for each category at that rank \
# If an ID list is specified in config write it   \
if config["method"] == "level":
    level = config["levelattr"]["level"].lower
    selectstr = '\'$rankId >= ' + str(leveldict[level]) + '\''
    reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree", "1", "--select", selectstr  ]
    reftreeout = get_reftree_data(reftreeopts)
    for line in reftreeout:
        trainlist, testlist = split_data(line, args.output, args.testfile,trainlist, testlist, float(config["testprop"]))

if config["method"] == "taxids":
    # Create a data list that contains taxa below selected nodes with the attributes specified in the config file
    sys.stderr.write("Formatting training data for selected taxa\n")
    # Search for records in each specified category
    for key, record in config["categories"].iteritems():  
        assert record["taxid"], "A taxon id is required for each category"
        
        # if gencode is specified, select only records with a matching gencode
        if "gencode" in record:
            selectstr = '\'$gencode == ' + str(record["gencode"]) + '\''    
            reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree", \
                str(record["taxid"]), "--select", selectstr]
            #print(reftreeopts)
        
        #otherwise return the the records without filtering 
        else:
            reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree", \
                str(record["taxid"])] 
            #print(reftreeopts)
        
        #Retrieve the records
        reftreeout = get_reftree_data(reftreeopts)
        # organelle is a read-level, not node-level attribute it must be parsed from the record after retrieval
        if reftreeout:
            if "organelle" in record:
                for line in reftreeout.split(os.linesep):
                    org = organelle(line)
                    if record["organelle"] == org:
                        trainlist, testlist = split_data(line, args.output, args.testfile,trainlist, testlist, float(config["testprop"]))  
                    else:
                        continue
            else:
                for line in reftreeout.split(os.linesep):
                    trainlist, testlist = split_data(line, args.output, args.testfile,trainlist, testlist, float(config["testprop"])) 
        
else:
    print("Unknown classification methods were requested in the config file.")
