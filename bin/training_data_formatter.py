#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to select taxonomic groups from a reftree feature vector database' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-o', '--training', help ="An output vector file", type=argparse.FileType('w'), default='-')
parser.add_argument('-v', '--crossvalidation', help ="An output vector file", type=argparse.FileType('w'), default = "crossval.txt")
parser.add_argument('-v', '--test', help ="An output vector file", type=argparse.FileType('w'), default = "crossval.txt")
parser.add_argument('-r', '--reftree', help ="Reftree database directory location") 
args = parser.parse_args()


# Functions
def get_reftree_data(dbpath, attributes):
    """A function to run RefTree on a particular taxonomy level or node and return the attributes and vector data"""
    reftreeopts = ["reftree.pl", "--dbDir", __ , "--db", __ , "--subtree", __ ]
    p0 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    reftreeout, reftreeerr= p0.communicate()
    assert p0.returncode == 0, "RefTree returned an error while searching the tree"
    return reftreeout

# Read the configuration file
config = json.load( open(args.config, 'r'))

test = float(config["split"]["test"])
train = float(config["split"["train"]) 
crossval = float(config["split"]["crossval"])
assert test + train + crossval = 1
	

# If an ID list is specified in config write it, otherwise query everything from the root node writing  \
if config["method"] == "level":
    print("Formatting training data for all taxa at the %s level" % config["levelattr"]["level"])
    level = config["levelattr"]["level"].lower
    data = get_reftree_data(config["dbpath"],level)
    args.output.write(data)
elif config["method"] == "taxids":
    print("Formatting training data for selected taxa")
    for key, record in config["categories"].iteritems():
        for pair in record:
            assert pair["taxid"], "A taxon id is required for each category"
            if "code" in pair:
                attributes = str(pair["taxid"]) + ":" + str(pair["code"])
            else:
                attributes = str(pair["taxid"])
            #print(str(key) + "\t" + attributes)
        	data = get_reftree_data(config["dbpath"],attributes)
        	args.output.write(data)
        	args.output.write("\n")
else:
    print("Unknown methods requested in hte config file.")
