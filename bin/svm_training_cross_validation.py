#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess


parser = argparse.ArgumentParser(description='A script for model training and validation from a Reftree database' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-r', '--reftree', help ="A local RefTree database of training vectors")
parser.add_argument('-o', '--output', help ="An output directory")
args = parser.parse_args()


tdfopts = ["training_data_formatter.py", "--config", args.config, "--reftree", args.reftree]
p0 = subprocess.Popen(tdfopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
tdfeout, tdferr= p0.communicate()
assert p0.returncode == 0, "RefTree returned an error while searching the tree"

