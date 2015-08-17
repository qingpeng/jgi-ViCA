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
def linetodict(line):
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

def mmatch(val):
    """returns a tuple containing the gencode and organelle from an item in the params dict"""
    if "organelle" in val:
        o = val["organelle"]
    else:
        o = None
    if "gencode" in val:
        g = val["gencode"]
    else:
        g = None
    return (o,g) 

def filter(d, params):
    # For each group in the config file
    for key, val in params.iteritems():
        taxid = val["taxid"] # this is the taxid of interest
        if taxid in d["taxonomy"]:
            a = mmatch(val)
            b = mmatch(d) 
            if a[0] and a[1]:
                if a[0] == b[0] and a[1] == b[1]:
                    return [key] + d["vector"]
            elif a[0]:
                if a[0] == b[0] and a[1] != b[1]:
                    return [key] + d["vector"]
            elif a[1]:
                if a[0] != b[0] and a[1] == b[1]:
                    return [key] + d["vector"]
            else:
                return [key] + d["vector"]
        else:
            continue


# MAIN

# retrieve all data
reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree" ,"1"] 
p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE)


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
    

