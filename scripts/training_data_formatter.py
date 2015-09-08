#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
import numpy


parser = argparse.ArgumentParser( description='A script to select taxonomic groups with particular attributes from a reftree feature vector database and split the data into test and train datasets' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-t', '--trainfile', help ="An output vector file", type=argparse.FileType('w'), default='-')
parser.add_argument('-e', '--testfile', help ="An output vector file", type=argparse.FileType('w'), default='testdata.txt')
parser.add_argument('-r', '--reftree', help ="Reftree database directory location") 
args = parser.parse_args()

# Read the configuration file
config = json.load(open(args.config, 'r'))["training_data_formatter"]

# Variables
dbdir,dbname = os.path.split(os.path.abspath(args.reftree))

# Functions

def validate(a,b):
    """A function to validate that every item in dictionary 'a' is present in 'b' and 
    its value matched the value in 'a' including cases where that value is 'None in 'a'
     but not present in 'b'. Returns True or False."""
    for key, val in a.iteritems():
        if val == "None":
            if key in b and b[key] != "None":
                return False
        else:
            if key not in b or b[key] != val:
                return False
    return True 
        
def linetodict(line):
    """split the data vector into a dictionary"""
    lv =line.strip().split("\t")
    if len(lv) < 5:
        return None
    c3 = lv[2].split(" ")
    attributes = lv[1]+","+ c3[2]
    attlist = attributes.split(',')
    d = dict(s.split('=') for s in attlist)
    taxlist = (d["taxonomy"].split('/'))
    t = dict(r.split(':') for r in taxlist)
    d["taxid"] = lv[0]
    d["taxonomy"] = t
    d["readId"] = c3[0]
    d["refseqId"] = c3[1]
    d["vector"] = lv[3:]
    return d


def filter(d, params):
    """Returns a class label and vector given  dictionaries of params and data"""
    # create a dictionary of lists containing the unique parameters for each taxid (needed to define wanted_keys)
    taxcatdict = {}
    for key1, val1 in params.iteritems():
        tmplist =[]
        for key2 in val1:
            tmplist.append(key2)
        tid = params[key1]["taxcat"]
        if tid in taxcatdict:
            taxcatdict[tid] = list(set(taxcatdict[tid] + tmplist))
        else:
            taxcatdict[tid] = tmplist 
    #iterate over the list pf parameters finding the category to assign the vector to   
    for key, val in params.iteritems():
        dc = d.copy()
        taxonomy = d["taxonomy"]
        if val["taxcat"] in taxonomy:
            dc["taxcat"] = val["taxcat"]
            if validate(val,dc) == True:
                dsub = {}
                wanted_keys = taxcatdict[val["taxcat"]] # The keys you want
                for item in wanted_keys:
                    if item in dc:
                        dsub[item] = dc[item]
                if validate(dsub, val) == True:
                    return [key] + dc["vector"]
                else:
                    continue
            else:
                continue
        else:
            continue
    

# MAIN

# validate config categories




# retrieve all data
reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree" ,"1"] 
p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE)

#


# Create lists to hold taxids
testlist = []
trainlist = []

# For each line 
for line in p0.stdout:
    # Convert the line to a dictionary
    d = linetodict(line)
    # decide whether data should go in the test or train bin
    if d:
        if d["taxid"] in trainlist:
            out = args.trainfile
        elif d["taxid"] in testlist:
            out = args.testfile
        else:
            selection = numpy.random.choice(["train","test"], p=[1-float(config["testprop"]),float(config["testprop"])])
            if selection == "train":
                out = args.trainfile
                trainlist.append(d["taxid"])
            elif selection == "test":
                out = args.testfile
                testlist.append(d["taxid"])
            else:
                raise ValueError('numpy.random.choice could not select a file to write the vector to')
        # Create a list with the category and vector if a match exists
        vectlist = filter(d, config["categories"])
        if vectlist: # If the vector list exists write it to the file
            out.write('\t'.join(vectlist))
            out.write('\n')
    else:
        continue
    

