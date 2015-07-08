#!/usr/bin/env python

import os
import argparse
import subprocess
import simplejson as json

parser = argparse.ArgumentParser(description='A script to extract genomic features on a single taxon \
	optionally shredding the taxon into sizes specified by a distribution and creating a \
	file with vector(s) describing the taxon and writing that taxon to a RefTree database')
parser.add_argument('-t', '--taxid', help ='An NCBI taxid') 
parser.add_argument('-i', '--input', help="A fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('-o', '--output', help ='A tab delimited list with taxid and a feature vector',\
	type=argparse.FileType('w'), default='-') 
parser.add_argument('-c', '--config', help='A JSON formatted configuration file', default='config.json')
args = parser.parse_args()

# Read the configuration file
config = json.loads( open(args.config, 'r'))

# set up scratch tmp dir
try:
	sge_task_id = os.environ['SGE_TASK_ID']
except:
	sge_task_id = 1
tempdir = os.path.join(config["nodetempdir"],'tmp_' + str(sge_task_id))
os.makedir(tempdir)

# Shred the sequence with shred.py
if conf['shred'] == 'gamma':
	shredopts = ["shred.py",  "--shred", "gamma","--samples",config["shredsamples"], \
		"--shape", config["gamma"]["shape"],"--scale", config["gamma"]["scale"], \
		"--offset", config["gamma"]["offset"]]
if conf['shred'] == 'fixed':
	shredopts = ["shred.py",  "--shred", "fixed", "--samples", config["shredsamples"], \
	"--length", config["gamma"]["length"]]
# If shredding is desired run shread.py
p1 = subprocess.Popen(shredopts, stdin=PIPE, stdout=PIPE)
p1.communicate(args.input) 
shredout, shrederr= p1.communicate()
	

#Run selected feature extraction script
if conf[method] == "metamark":
	#run metamark wrapper
	metamarkopts = ["feature_extraction_metamark.py", "-config", config, "-tmpdir", tempdir]
	p2 = subprocess.Popen(metamarkopts, stdin=p1.PIPE, stdout=PIPE)
	p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
	matrixdata, metamarkwerr= p2.communicate()
if conf[method] == "kmer":
	print("kmer method not yet implemented")
	pass

# save vector to reftree
for line in matrixdata:
	args.output.write(line + "\n")
 	

 	
