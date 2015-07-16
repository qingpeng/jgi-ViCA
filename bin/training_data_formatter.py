#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script select taxonomic groups from a reftree feature vector database' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-o', '--output', help ="An output vector file", )
parser.add_argument('-a', '--attribute', help ="Reftree database directory location") 
args = parser.parse_args()


# Functions

def get_reftree_data(dbpath, attributes):
	"""A function to run reftree on a particular taxonomy level and return the class and vector data"""
	reftreeopts = ["reftree.pl", "--tree", "dbpath", attributes]
	p0 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	reftreeout, reftreeerr= p0.communicate()
	assert p0.returncode == 0, "RefTree returned an error while searching the tree"
	return reftreeout



# Read the configuration file
config = json.load( open(args.config, 'r'))


# If an id list is specified in config write it otherwise query everything from the root node \
if config["method"] == "level":
	print("Formatting training data for all taxa at the %s level" % config["levelattr"]["level"]
	level = config["levelattr"]["level"].lower
	data = get_reftree_data(config[dbpath],level)
elif config["method"] == "taxids":
	print("Formatting training data for selected taxa")
	for key, record in Config["Categories"]:
	TODO
else:
	print("Unknown methods requested in hte config file.")
