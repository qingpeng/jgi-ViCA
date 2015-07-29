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


# Read the configuration file
config = json.load(open(args.config, 'r'))["training_data_formatter"]

dbdir,dbname = os.path.split(os.path.abspath(args.reftree))

# Functions

#todo add level and attribute filtering once Ed adds it to Reftree
def get_reftree_data(reftreeopts):
    """A function to run RefTree"""
    p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE)
    reftreeout, reftreeerr= p0.communicate()
    assert p0.returncode == 0, "RefTree returned an error while searching the tree with taxid"
    return reftreeout

        


# If a a taxonomic rank is specified return the data for for each category at that rank \
# If an ID list is specified in config write it   \
# if config["method"] == "level":
#   pass
#     sys.stderr.write("Formatting training data for all taxa at the %s level" % config["levelattr"]["level"])
#     level = config["levelattr"]["level"].lower
#     reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--some_command?", level ]
#     data = get_reftree_data(reftreeopts)
#     args.output.write(data)
#     #To be done after Ed implements level feature
if config["method"] == "taxids":
    # Create a data list that contains taxa below selected nodes with the attributes specified in the config file
    sys.stderr.write("Formatting training data for selected taxa\n")
    for key, record in config["categories"].iteritems():  
        assert record["taxid"], "A taxon id is required for each category"
        attributes = [str(record["taxid"]), "--attributes", "code"]
        reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree"] + attributes
        print(reftreeopts)
#         data = get_reftree_data(reftreeopts)
#         data_list = data.split()
#         if record["code"]:
#             if data_list[1] == record["code"]:
#                 final = [key] + data_list[2:]
#         else:
#             final = [key] + data_list[2:]
#         args.output.write(final)
#         args.output.write("\n")
else:
    print("Unknown classification methods were requested in the config file.")
