#!/usr/bin/env python
#create_vectors.py
#adam Rivers

import argparse
import os
import subprocess
#import sys
from itertools import chain
import cPickle as pickle

parser = argparse.ArgumentParser(description='A script to generate SVM training matrices from Genemark .mod files')
parser.add_argument('-r','--root', help='a root input directory', required=True)
parser.add_argument('-c,','--config',help= 'a json formated config file', default='config.json')
args = parser.parse_args()


root = os.path.abspath(args.root)	
mod = os.path.join(root, "data","trainingmod")

def flatten(listOfLists):
    """Flatten one level of nesting"""
    return chain.from_iterable(listOfLists)

classvect = []
modvect = []

for path, dirs, files in os.walk(mod, topdown=False):
	for filename in files:
		fullpath = os.path.join(path, filename)
		if os.path.isfile(fullpath):
			with open(fullpath, 'r') as f:
				if f.name.endswith("hmm.mod"):
 					rawvec = []
 					for line in f.readlines():
						rawvec.append(line.split())
					start1 = rawvec.index(['COD1']) +1
					stop1 = rawvec.index(['COD2']) -1
					COD1 = (list(flatten(rawvec[start1: stop1 ])))
					start2 = 1 + rawvec.index(['NONC'])
					NONC = (list(flatten(rawvec[start2: start2 +64 ])))
					itemvect = COD1 + NONC
					itemvect = [float(i) for i in itemvect]
					assert len(itemvect) == 256, "a full length vector could not be extracted from %s" % f.name
					modvect.append(itemvect)
					classvect.append(os.path.basename(path))
					 		
#print(modvect)
#print(classvect)
outfile = os.path.join(root,"SVMmatrix.p")
with open(outfile, 'wb' ) as outf:
	pickle.dump(modvect, outf)
	pickle.dump(classvect, outf)

print("completed run")