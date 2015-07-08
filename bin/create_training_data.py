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
parser.add_argument('-o', '--output', help ="Reftree database directory") 
args = parser.parse_args()


# Read the configuration file
config = json.loads( open(args.config, 'r'))
# If an id list is specified in config write it otherwise query everything from the root node \
if config["trainingid"]:
	nodes = config["trainingid"]
	print("Loading Subtree taxon ids from config file: %i" % nodes)
else:
	nodes = [1]
	print("No subtree specified, training on all sequences")
	
# Create data directory if necessary and 100 subdirectories
datadir = os.path.join(os.path.abspath(config["rootdir"],"data"))
if not os.path.exists(datadir):
	os.makedirs(datadir)
for i in range(10,100):
	os.makedirs(os.path.join(datadir,str(i))

# Open a subprocess to have RefTree pass taxids
taxids = []
for node in nodes:
	reftreeparams = ["refrree.pl","--taxonomy", node, "genomic" ]
	p1 = subprocess.Popen(reftreeparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	reftreeout, reftreeerr= p1.communicate()
	taxids.append(reftreeout)
	assert p1.returncode == 0, "RefTree returned the error %s while looking up the \
		taxonomy ids" % qsuberr
	# root/data/01/0143242.fasta, etc

## Create generate a task farmer job list file ##

# Create a Task farmer directory
tfdir = os.path.join(os.path.abspath(config["analysisrootdir"]),"taskfarmer")
if not os.path.exists(tfdir):
	os.path.makedir(tfdir)

# Open the Task Farmer list file
tflist = os.path.join(ftdir, "taskfarmer.list")
with open(tflist, 'w') as tfl: 
	# for every fasta write a line in the task farmer list file
	for taxid in taxids:
		outfile = os.path.abspath(taxid) +".done"
		tfl.write("single_taxon_training.py "--taxid"+ taxid + ":" +  outfile + ":" + "0")
tflist.close()

# Start up the task farmer "client" i.e. manager
p2 = subprocess.Popen(["tfmq-client", "-i", "tflist" "-q", "singletaxontraining"], \
	stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
tfmqout,tfmqerr = p2.communicate()

# qsub the script containing the Task Farmer workers
p3 = subprocess.Popen(["qsub", "-l", "ram.c=5G", "-t" "1-40", "./taskfermer_sts.sh"], \
	stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
qsubout,qsuberr = p2.communicate()
assert p3.returncode == 0, "qsub returned the error %s" % qsuberr

# Monitor the Task farmer client and wait for the processes to complete
if p2.returncode == 0:
	
	# Call Reftree to have it generate reftree local db from the datadir
	reftreedb = os.path.abspath(args.output)
	if not os.path.exists(reftreedb):
		os.makedb(reftreedb)
	p4 = subprocess.Popen(["reftree.pl", "--loadDb", datadir, reftreedb ], \
		stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	tr2out,tr2err = p2.communicate()