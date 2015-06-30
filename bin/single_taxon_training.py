#!/usr/common/usg/languages/python/2.7.4/bin/python

import os
import argparse
import subprocess
import simplejson as json
import jsonschema
import random
import string
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to run genmark on a single taxon \
	optionally shredding the taxon into sizes spefified by a distbution and creating a \
	file with vector(s) describing the taxon')
parser.add_argument('-i', '--input', type=argparse.FileType('r'), help="A fasta file",\
	 default='-')
parser.add_argument('-t', '--taxid', help ="An NCBI taxid") 
parser.add_argument('-c','--config', help="A JSON formatted configuration file")
parser.add_argument('-s','--schema', help="A JSON formatted configuration schema", \
	default="config_schema.json")
args = parser.parse_args()


#Parse config file
conf = json.load(args.config)
#validation to be added latter
# confschema = json.load(args.schema)
# try:
# 	validate(conf, confschema)
# except ValidationError:
# 	print("There seems to be a problem with the format of your configuration file")



#get root directory
root = os.path.abspath(config["rootdir"])
temp ='temp_' + ''.join(random.choice(string.ascii_uppercase) for _ in range(6))
taxonfile = os.path.join(root,args.taxid)
if not os.path.exists(taxondir):
    os.makedirs(d=taxondir)
tempdir = os.path.join(taxondir,temp)
if not os.path.exists(tempdir)
	os.makedir(tempdir)



if config['shread'] == gamma:
	grinderopts = ["grinder", "-rf", "-" "-cf" "0.1", "-rd", config["gamma"]["a"] \
		,"gamma", config["gamma"]["b"], "-un", "-od", tempdir]
if config['shread'] = fixed:
	grinderopts = ["grinder", "-rf", "-" "-cf" "0.1", "-rd", config["fixed"], "uniform"\
		 ,"0", "-un", "-od", tempdir]
	#If shredding is desired run grinder
	p1 = subprocess.Popen(grinderopts, stdin=PIPE, stdout=subprocess.PIPE)
		grinderout, grindererr= p1.communicate()
		p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.

#Run selected feature extraction script
if conf[method] == "metamark":
	pass
if conf[method] == "kmer":
	pass

# convert csv data to json with matrix_2_json.py

# add to refTree


