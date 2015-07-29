#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
import sys


parser = argparse.ArgumentParser(description='A script to select taxonomic groups from a reftree feature vector database' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-o', '--output', help ="An output vector file", type=argparse.FileType('w'), default='-')
parser.add_argument('-r', '--reftree', help ="Reftree database directory location") 
args = parser.parse_args()


# Functions

#todo add level and attribute filtering once Ed adds it to Reftree
def get_reftree_data(dbdir, dbname, subtree, level, code):
    """A function to run RefTree on a particular taxonomy level or node and return the attributes and vector data"""
    if level:
    	reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname ,  level ]
    else:
    	reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree", subtree ]
    p0 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    reftreeout, reftreeerr= p0.communicate()
    assert p0.returncode == 0, "RefTree returned an error while searching the tree with taxid %s" % subtree
    return reftreeout

# Read the configuration file
config = json.load( open(args.config, 'r'))


# If a a taxonomic rank is specified return the data for for each category at that rank \
# If an ID list is specified in config write it   \
if config["method"] == "level":
    sys.stderr.write("Formatting training data for all taxa at the %s level" % config["levelattr"]["level"])
    level = config["levelattr"]["level"].lower
    data = get_reftree_data(dbdir = config["dbpath"],level = level)
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
    print("Unknown classification methods were requested in the config file.")
