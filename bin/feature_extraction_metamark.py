#!/usr/bin/env python

import sys
import simplejson as json
import argparse
import os
import subprocess
from Bio import SeqIO
from itertools import chain
import csv
import shutil

parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
using emmission data from Metamark')
parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help= "output file csv format", type=argparse.FileType('w'), default='-')
parser.add_argument('--tmpdir', help="A path to a temporary directory")
parser.add_argument('--taxid', help="A path to a temporary directory")
parser.add_argument('--mmp', help="the parameters file for metamark", default = "../gm_parameters/par_11.modified")
parser.add_argument('-c','--config', help="A JSON formatted configuration file")

args = parser.parse_args()

## File parsing and variable assignment
path = os.path.join(os.path.abspath(args.tmpdir),"metamark" )
records = SeqIO.parse(args.input, "fasta")
mmp = os.path.abspath(args.mmp)
csvout = csv.writer(args.output, quoting=csv.QUOTE_MINIMAL)


# Functions
def flatten(listOfLists):
    """Flatten one level of nesting"""
    return chain.from_iterable(listOfLists)

def parsemod(dir, taxid):
    # Find file and parse
    for file in os.listdir(dir):
        if file.endswith("hmm.mod"):
            with open(file, "r") as f:
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
                itemvect[:0] = str(taxid)
                ### This isn't working: splits taxid into digits
                print(taxid)
                return itemvect


# for each sequence in the fasta file:
for record in records:
    #os.mkdir(path) # Generate a temp dir
    os.chdir(path) 
    handle = open("fragment.fasta", "w") # open a fasta file
    SeqIO.write(record, handle, "fasta") # write the sequence to it
    handle.close() # close the file
    
    ## Run metamark
    metamarkparams = ["gmsn.pl", "--clean", "--gm", "--par", mmp,"fragment.fasta"]
    p1 = subprocess.Popen(metamarkparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    metamarkout, metamarkerr= p1.communicate()
    #print(metamarkout)
    #print(metamarkerr)
    #Parse results
    dir = os.getcwd()
    featurevect = parsemod(dir,args.taxid)
    print(featurevect)
#     csvout.writerow(featurevect)
#     os.chdir(os.path.abspath(args.tmpdir))
#     shutil.rmtree(path)
    
