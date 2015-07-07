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

#functions
def parsegenomes(inhandle, root):
	"""Splits a stream of genomes from Reftree into individual genome files"""
	sequences = SeqIO.parse(inhandle, 'fasta')
	for sequence in sequences:
		taxid = re.sub('.\|tax\|','',sequence.id)
		geneid = taxid = re.sub('\|tax\|.','',sequence.id)
		genomedir = os.path.join(root,str(taxid)[0:2])
		if not os.path.exists(genomedir):
			os.makedirs(genomedir)
		filename = re.sub('\|','\-',sequence.id) + "fasta"
		with open(os.path.join(genomedir,filename), 'w') as outhandle:
			SeqIO.write(sequence,outhandle,"fasta")
	outhandle.close()



# Read config file
config = json.loads( open(args.config, 'r'))
# If an id list is specified in config write it otherwise query everything from the root node \
if config["trainingid"]:
	nodes = config["trainingid"]
	print("Loading Subtree taxon ids from config file: %i" % nodes)
else:
	nodes = [1]
	print("No subtree specified, training on all sequences")
	
#query reftree for genomes writing them to a directory specified in config
# create data directory if necessary
datadir = os.path.join(os.path.abspath(config["analysisrootdir"],"data"))
if not os.path.exists(datadir):
	os.makedirs(datadir)
# open a subprocess to have reftree pass genome files
for node in nodes:
	reftreeparams = ["RefTree.pl","--tree", node, "genomic" ]
	p1 = subprocess.Popen(reftreeparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	reftreeout, reftreeerr= p1.communicate()
	parsegenomes(reftreeout,datadir)
	# root/data/01/0143242.fasta, etc

# Create generate a task farmer job list file
#create a Task farmer directory
tfdir = os.path.join(os.path.abspath(config["analysisrootdir"]),"taskfarmer")
if not os.path.exists(tfdir):
	os.path.makedir(tfdir)
# open the list file
tflist = os.path.join(ftdir, "taskfarmer.list")
with open(tflist, 'w') as tfl: 
	# for every fasta write a line in the task farmer list file
	for file in os.path.walk(datadir):
		taxid = re.sub('.-tax-','',sequence.id)
		outfile = os.path.abspath(file) +".done"
		tfl.write("single_taxon_training.py --input " + file + "--taxid"+ taxid + ":" +  outfile + ":" + "0")
tflist.close()

#initiate the task farmer client
p2 = subprocess.Popen(["tfmq-client", "-i", "tflist" "-q", "singletaxontraining"], \
	stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
tfmqout,tfmqerr = p2.communicate()

# qsub the workers
p3 = subprocess.Popen(["qsub", "-l", "ram.c=5G", "-t" "1-40", "./taskfermer_sts.sh"], \
	stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
qsubout,qsuberr = p2.communicate()
assert p3.returncode == 0, "qsub returned the error %s" % qsuberr

# monitor the Task farmer client and wait for the processes to complete
if p2.returncode == 0:
	#call Reftree to have it generate reftree local db from the datadir
	reftreedb = os.path.abspath(args.output)
	if not os.path.exists(reftreedb):
		os.makedb(reftreedb)
	p4 = subprocess.Popen(["reftree.pl", "--loadDb", datadir, reftreedb ], \
		stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	tr2out,tr2err = p2.communicate()