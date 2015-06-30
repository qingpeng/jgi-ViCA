#!/usr/common/usg/languages/python/2.7.4/bin/python

import os
import argparse
import subprocess
import simplejson as json
import jsonschema
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to run genmark on a single taxon optionally shredding the taxon into sizes spefified by a distbution and creating a file with vector(s) describing the taxon')
parser.add_argument('-i', '--input', type=argparse.FileType('r'), help="A fasta file", , default='-')
parser.add_argument('-c','--config', help="a JSON formatted configuration file")
parser.add_argument('-s','--schema', help="a JSON formatted configuration file", default="config_schema.json")
args = parser.parse_args()

#Parse config file
conf = json.load(args.config)
confschema = json.load(args.schema)
try:
	validate(conf, confschema)
except ValidationError:
	print("There seems to be a problem with the format of your configuration file")

#If shredding is desired run Grinder

#Run selected feature extraction script
if conf[method] == "metamark":
	pass
if conf[method] == "kmer":
	pass

# convert csv data to json with matrix_2_json.py

# add to refTree


