#!/usr/bin/env python

import os
import argparse
import subprocess
import simplejson as json
import shutil
import sys
import signal


parser = argparse.ArgumentParser(description='A script to extract genomic features on a single taxon \
    optionally shredding the taxon into sizes specified by a distribution and creating a \
    file with vector(s) describing the taxon and writing that taxon to a RefTree database')
parser.add_argument('-t', '--taxid', help ='An NCBI taxid',required=True)
parser.add_argument('-o', '--outfile', help ='output file with shredded segments', required=True)
parser.add_argument('-c', '--config', help='A JSON formatted configuration file', default='config.json')
args = parser.parse_args()

def signal_term_handler(signal, frame):
	sys.stderr.write('single_taxon_training.py was killed by a SIGTERM command')
	sys.exit(1)

signal.signal(signal.SIGTERM, signal_term_handler)

# Read the configuration file
config = json.load( open(args.config, 'r'))["create_training_data"]

## Retrieve the sequence from reftree
# with taxonomy -QP for reading raw results
reftreeopts = ["reftree.pl" , "--db", "genomic", "--node", args.taxid,"--attributes", "gencode,taxonomy"]
p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE, stderr = subprocess.PIPE)

output_file = args.outfile


## Shred the sequence with shred.py
if config['shred'] == 'lognorm':
    shredopts = ["shred.py",  "--shred", "lognorm","--samples",config["shredsamples"], \
        "--shape", config["lognorm"]["shape"],"--scale", config["lognorm"]["scale"], \
        "--loc", config["lognorm"]["loc"], "--output", output_file]
elif config['shred'] == 'fixed':
    shredopts = ["shred.py",  "--shred", "fixed", "--samples", config["shredsamples"], \
    "--length", config["fixed"], "--output", output_file]
#    shredopts2 = ["shred.py",  "--shred", "fixed", "--samples", "200", \
#    "--length", "5000", "--output", "shred2.fa","--input","138_genome.fa"]
#    shredopts3 = ["shred.py",  "--shred", "fixed", "--samples", "200", \
#    "--length", "5000", "--output", "shred3.fa"]
    
    
elif config['shred'] =='None':
	output_file_obj = open(output_file,'w')
	output_file_obj.write(p0.stdout)
	
else:
	raise ValueError('An inccorect value was specified in the configuration file for shredding: ' + config['shred'])

if config['shred'] == 'lognorm' or config['shred'] == 'fixed':
	# If shredding is desired run shread.py
    print "fixed"
    p1 = subprocess.Popen(shredopts, stdin=p0.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#p2 = subprocess.Popen(shredopts2, stdin=p0.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#for line in  p0.stdout:
	#    print line
	    
#    p3 = subprocess.Popen(shredopts3, stdin=p0.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#p3 = subprocess.Popen(shredopts3, stdin=p0.stdout)
	#for line in p3.stdout:
	#    print line
#	p3.stdout.close()
    matrixdata, metamarkerr= p1.communicate()

    p0.stdout.close()  #This is needed in order for p0 to receive a SIGPIPE if p1 exits before p0

