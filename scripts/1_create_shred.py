#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to create sequences for training from gemomes' )
parser.add_argument('-c', '--config', help ="json formatted config file", required= True)
parser.add_argument('-o', '--output', help ="directory to store the shredded sequences",required= True) 
args = parser.parse_args()


# Variables
configpath = os.path.abspath(args.config)

# Read the configuration file
config = json.load(open(configpath, 'r'))["create_training_data"]
config_sys = json.load(open(configpath, 'r'))["environment_variables"]

sts = os.path.join(str(config_sys["genelearn_path"]),"scripts","single_taxon_shred.sh")


# If an id list is specified in config write it otherwise query everything from the root node \
if config["trainingid"]:
	node = str(config["trainingid"])
	print("Loading subtree taxon id from config file: %s" % node)
else:
	node = str(1)
	print("No subtree specified, training on all sequences")
	
# Open a subprocess to create a new reftree database directory, have it execute a shell \
# script that calls single_taxon_training and writes the results to the reftree database directory


resources = config["resources"]
reftreeopts = ["reftree.pl", "--db", "genomic","--resources", resources, "--keep", "--foreach", node, args.output,\
"--",sts, "__TAXON__","__OUTFILE__", configpath]
#print(reftreeopts)
p1 = subprocess.Popen(reftreeopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
reftreeout, reftreeerr= p1.communicate()
#print(reftreeout)

