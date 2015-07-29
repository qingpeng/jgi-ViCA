#!/usr/bin/env python

import os
import argparse
import subprocess
import simplejson as json
import biopython


parser = argparse.ArgumentParser(description='A script to generate an ESOM data file from contigs')
parser.add_argument('-t', '--infile', help ='a file with contigs',type=argparse.FileType('r'), default='-') 
parser.add_argument('-o', '--outfile', help ='location to write the feature vector file', type=argparse.FileType('w'), default='-') 
parser.add_argument('-c', '--config', help='A JSON formatted configuration file', default='config.json')
args = parser.parse_args()

config = json.load(open(configpath, 'r'))["create_training_data"]

femopts = ["feature_extraction_metamark.py","--label", "readid", "--mmp", config["mmp"] ]
p1 = subprocess.Popen(femopts, stdin=args.infile, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#assert p1.returncode == 0, "feature_extraction_metamark.py encountered an error"


ffopts = ['feature_formatter.py', '--infmt','csv','-outfmt','esom-lrn']
# If shredding is desired run shread.py
p1 = subprocess.Popen(ffopts, stdin=p0.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p0.stdout.close()  #This is needed in order for p0 to receive a SIGPIPE if p1 exits before p0
ffout, fferr = p1.communicate()
assert p1.returncode ==0, "the feature formatter script could not return a valid file"
args.outfile.write(ffout)