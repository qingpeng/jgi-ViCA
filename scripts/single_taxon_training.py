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
parser.add_argument('-o', '--outfile', help ='location to write the  feature vector file', required=True) 
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



## Shred the sequence with shred.py
if config['shred'] == 'lognorm':
    shredopts = ["shred.py",  "--shred", "lognorm","--samples",config["shredsamples"], \
        "--shape", config["lognorm"]["shape"],"--scale", config["lognorm"]["scale"], \
        "--loc", config["lognorm"]["loc"]]
elif config['shred'] == 'fixed':
    shredopts = ["shred.py",  "--shred", "fixed", "--samples", config["shredsamples"], \
    "--length", config["fixed"]]
elif config['shred'] =='None':
	sequencein = p0.stdout 
else:
	raise ValueError('An inccorect value was specified in the configuration file for shredding: ' + config['shred'])

if config['shred'] == 'lognorm' or config['shred'] == 'fixed':
	# If shredding is desired run shread.py
	print "fixed"
	p1 = subprocess.Popen(shredopts, stdin=p0.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	p0.stdout.close()  #This is needed in order for p0 to receive a SIGPIPE if p1 exits before p0
	sequencein = p1.stdout 
#temp debug
print "test"
#for line in sequencein:
#    print line

print sequencein
## Run selected feature extraction script
if config["method"] == "genemarks":
    #run metamark wrapper
    genemarksopts = ["feature_extraction_genemarks.py", "--tmp", config["tmpdir"],"--mmp", config["mmp"], "--taxid", args.taxid, "--outfile", args.outfile]
    p2 = subprocess.Popen(genemarksopts, stdin=sequencein , stdout=subprocess.PIPE)
    p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
    matrixdata, metamarkerr= p2.communicate()
    assert p2.returncode == 0, 'there was an error in single taxon training with taxid %s' % args.taxid

elif config["method"] == "metagenemark":
    #run metamark wrapper
    metagenemarkopts = ["feature_extraction_metagenemark.py", "--tmp", config["tmpdir"],"--mmp", config["mmp"], "--taxid", args.taxid, "--outfile", args.outfile]
    p2 = subprocess.Popen(metagenemarkopts, stdin=sequencein , stdout=subprocess.PIPE)
    p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
    matrixdata, metamarkerr= p2.communicate()
    assert p2.returncode == 0, 'there was an error in single taxon training with taxid %s' % args.taxid

elif config["method"] == "metagenemark_kmer":
    #run metamark wrapper
    metagenemarkopts = ["feature_extraction_metagenemark_kmer.py", "--tmp", config["tmpdir"],"--mmp", config["mmp"], "--taxid", args.taxid, "--outfile", args.outfile]
    p2 = subprocess.Popen(metagenemarkopts, stdin=sequencein , stdout=subprocess.PIPE)
    p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
    matrixdata, metamarkerr= p2.communicate()
    assert p2.returncode == 0, 'there was an error in single taxon training with taxid %s' % args.taxid
    
elif config["method"] == "kmer":
    metamarkopts = ["feature_extraction_kmer.py", "--taxid", args.taxid, "--outfile", args.outfile]
    p2 = subprocess.Popen(metamarkopts, stdin=sequencein , stdout=subprocess.PIPE)
    p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
    matrixdata, metamarkerr= p2.communicate()
    assert p2.returncode == 0, 'there was an error in single taxon training with taxid %s' % args.taxid

elif config["method"] == "both":
    #run metamark wrapper
    metamarkopts = ["feature_extraction.py", "--tmp", config["tmpdir"],"--mmp", config["mmp"], "--taxid", args.taxid, "--outfile", args.outfile]
    p2 = subprocess.Popen(metamarkopts, stdin=sequencein , stdout=subprocess.PIPE)
    if config['shred'] != 'None':
        p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
    matrixdata, metamarkerr= p2.communicate()
    assert p2.returncode == 0, 'there was an error in single taxon training with taxid %s' % args.taxid

elif config["method"] == "metagenemark_v2":
    #run metamark wrapper
    metamarkopts = ["feature_extraction_v2.py", "--tmp", config["tmpdir"],"--mmp", config["mmp"], "--taxid", args.taxid, "--outfile", args.outfile]
    p2 = subprocess.Popen(metamarkopts, stdin=sequencein , stdout=subprocess.PIPE)
    if config['shred'] != 'None':
        p1.stdout.close()  #This is needed in order for p1 to receive a SIGPIPE if p2 exits before p1
    matrixdata, metamarkerr= p2.communicate()
    assert p2.returncode == 0, 'there was an error in single taxon training with taxid %s' % args.taxid

