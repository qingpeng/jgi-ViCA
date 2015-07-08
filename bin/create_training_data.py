#!/usr/bin/env python

import argparse
import simplejson as json
import jsonschema
import csv
import os
import subprocess
from Bio import SeqIO
parser = argparse.ArgumentParser(description='A script to create training data from gemomes' )
parser.add_argument('-c', '--config', help ="json formatted config file")
parser.add_argument('-o', '--output', help ="Reftree database directory location") 
args = parser.parse_args()


# Read the configuration file
config = json.loads( open(args.config, 'r'))
# If an id list is specified in config write it otherwise query everything from the root node \
if config["trainingid"]:
	node = config["trainingid"]
	print("Loading subtree taxon ids from config file: %i" % nodes)
else:
	node = [1]
	print("No subtree specified, training on all sequences")
	
# Open a subprocess to create a new reftree database directory, have it execute a shell \
# script that calls single_taxon_training and writes the results to the reftree database directory

reftreeopts = ["reftree.pl", "--newdb", args.output, "--submit", "genomic", \
	"./single_taxon_training.py","TRUE","--","__TAXON__","__TMPFILE__"]
p1 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
reftreeout, reftreeerr= p1.communicate()
assert p1.returncode == 0, "RefTree returned the error %s while executing taxon training" % qsuberr


# Index the reftree database #
reftree2opts = ["reftree.pl", "--indexDb", args.output]
p2 = subprocess.Popen(reftree2opts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
reftreeout, reftreeerr= p1.communicate()
assert p2.returncode == 0, "RefTree returned the error %s while indexing the database" % qsuberr
