#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to create training data from gemomes' )
parser.add_argument('-c', '--config', help ="json formatted config file")
parser.add_argument('-o', '--output', help ="Reftree database directory location") 
args = parser.parse_args()


#Variables
configpath = os.path.abspath(args.config)
# Read the configuration file
config = json.load(open(configpath, 'r'))
numworkers = str(config["numworkers"])
resources = config["resources"]
# If an id list is specified in config write it otherwise query everything from the root node \
if config["trainingid"]:
	node = str(config["trainingid"])
	print("Loading subtree taxon id from config file: %s" % node)
else:
	node = str(1)
	print("No subtree specified, training on all sequences")
	
# Open a subprocess to create a new reftree database directory, have it execute a shell \
# script that calls single_taxon_training and writes the results to the reftree database directory

reftreeopts = ["reftree.pl", "--db", "genomic","--slots",numworkers, "--resources", resources,"--keep", "--foreach", node, args.output,\
"--","/global/homes/a/arrivers/dev/jgi-genelearn/bin/single_taxon_training.sh", "__TAXON__","__OUTFILE__", configpath]
#print(reftreeopts)
p1 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
reftreeout, reftreeerr= p1.communicate()
print(reftreeout)

