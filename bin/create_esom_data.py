#!/usr/bin/env python

import os
import argparse
import subprocess
import simplejson as json
import sys

parser = argparse.ArgumentParser(description='A script to generate an ESOM data file from contigs')
parser.add_argument('-', '--infile', help ='a file with contigs',required=True) 
parser.add_argument('-o', '--outfile', help ='location to write the feature vector file', type=argparse.FileType('w'), default='-') 
parser.add_argument('-c', '--config', help='A JSON formatted configuration file', default='config.json')
args = parser.parse_args()

config = json.load(open(args.config, 'r'))["create_training_data"]
infile = os.path.abspath(args.infile)
femopts = ["feature_extraction_metamark.py", "--input", infile, "--label", "readid", "--mmp", config["mmp"] ]
p1 = subprocess.Popen(femopts, stdout=subprocess.PIPE,stderr=subprocess.PIPE)

ffopts = ['feature_formatter.py', '--infmt','csv','--outfmt','esom-lrn']
p2 = subprocess.Popen(ffopts, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p1.stdout.close()  #This is needed in order for p0 to receive a SIGPIPE if p1 exits before p0
ffout, fferr = p2.communicate()
p1.wait()
assert p2.returncode ==0, "the feature formatter script could not return a valid file"
sys.stderr.write(fferr)
args.outfile.write(ffout)
