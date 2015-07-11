#!/usr/bin/env python

import os
import argparse
import subprocess
import simplejson as json
import shutil

parser = argparse.ArgumentParser(description='A script to extract genomic features on a single taxon \
	optionally shredding the taxon into sizes specified by a distribution and creating a \
	file with vector(s) describing the taxon and writing that taxon to a RefTree database')
parser.add_argument('-t', '--taxid', help ='An NCBI taxid',required=True) 
parser.add_argument('-i', '--input', help="A fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('-o', '--outfile', help ='location to write the  feature vector file', required=True) 
parser.add_argument('-c', '--config', help='A JSON formatted configuration file', default='config.json')

args = parser.parse_args()

# Read the configuration file
config = json.load( open(args.config, 'r'))

# set up scratch tmp dir
try:
	sge_task_id = os.environ['SGE_TASK_ID']
except:
	sge_task_id = 1
tempdir = os.path.join(config["nodetempdir"],'tmp_' + str(sge_task_id))
os.makedirs(tempdir)

## Shred the sequence with shred.py
if config['shred'] == 'gamma':
	shredopts = ["shred.py",  "--shred", "gamma","--samples",config["shredsamples"], \
		"--alpha", config["gamma"]["alpha"],"--scale", config["gamma"]["scale"], \
		"--loc", config["gamma"]["loc"]]

if config['shred'] == 'fixed':
	shredopts = ["shred.py",  "--shred", "fixed", "--samples", config["shredsamples"], \
	"--length", config["fixed"]]

# If shredding is desired run shread.py
p1 = subprocess.Popen(shredopts, stdin=args.input, stdout=subprocess.PIPE)

#Run selected feature extraction script
if config["method"] == "metamark":
	#run metamark wrapper
	metamarkopts = ["feature_extraction_metamark.py", "--tmpdir", tempdir, "--mmp", config["mmp"], "--taxid", args.taxid, "--outfile", args.outfile]
	p2 = subprocess.Popen(metamarkopts, stdin=p1.stdout , stdout=subprocess.PIPE)
	p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
	matrixdata, metamarkerr= p2.communicate()

elif config["method"] == "kmer":
 	print("kmer method not yet implemented")
 	
shutil.rmtree(tempdir)