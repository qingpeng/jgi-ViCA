#!/usr/bin/env python

import os
import argparse
import subprocess
import simplejson as json
import sys
import tempfile
import shutil

parser = argparse.ArgumentParser(description='A script to generate an ESOM data file from contigs')
parser.add_argument('-i', '--infile', help ='a file with contigs',required=True) 
parser.add_argument('-o', '--outfile', help ='location to write the feature vector file', required = True) 
parser.add_argument('-c', '--config', help='A JSON formatted configuration file', default='config.json')
args = parser.parse_args()

tmpdir = tempfile.mkdtemp(dir="/scratch")
tmpfile = os.path.join(tmpdir, "vector.txt")

config = json.load(open(args.config, 'r'))["create_training_data"]
infile = os.path.abspath(args.infile)
femopts = ["feature_extraction_metamark.py", "--input", infile, "--outfile", tmpfile, "--label", "readid", "--mmp", config["mmp"] ]
p1 = subprocess.Popen(femopts, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
femout,femerr = p1.communicate()
if p1.returncode != 0:
	shutil.rmtree(tmpdir)
	print(femout)
	print(femerr)
	raise RunError("Feature extraction failed")

ffopts = ['feature_formatter.py', "--input",tmpfile , "--output",args.outfile, '--infmt','csv','--outfmt','esom-lrn']
p2 = subprocess.Popen(ffopts, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
p1.stdout.close()  #This is needed in order for p0 to receive a SIGPIPE if p1 exits before p0
ffout, fferr = p2.communicate()
if  p2.returncode != 0:
	shutil.rmtree(tmpdir)
	print(ffout)
	print(fferr)
	raise RunError("the feature formatter script could not return a valid file")

shutil.rmtree(tmpdir)