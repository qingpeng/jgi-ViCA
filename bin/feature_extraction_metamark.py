#!/usr/bin/env python
import argparse
import os
import subprocess
from Bio import SeqIO
from itertools import chain
import shutil
import re

parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
using emmission data from Metamark')
parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
parser.add_argument('--tmpdir', help="A path to a temporary directory")
parser.add_argument('--taxid', help="The taxonomy id")
parser.add_argument('--mmp', help="the parameters file for metamark", default = "../gm_parameters/par_11.modified")
args = parser.parse_args()

## File parsing and variable assignment
path = os.path.join(os.path.abspath(args.tmpdir),"metamark" )
records = SeqIO.parse(args.input, "fasta")
mmp = os.path.abspath(args.mmp)

# Functions
def flatten(listOfLists):
    """Flatten a list to one level of nesting"""
    return chain.from_iterable(listOfLists)

def parsemod(dir, taxid):
    """Finds the Genemark model, parses it and returns a feature vector with the taxid as the first element"""
    for file in os.listdir(dir):
        if file.endswith("hmm.mod"):
            with open(file, "r") as f:
                rawvec = []
                taxlist = [taxid]
                for line in f.readlines():
                    rawvec.append(line.split())
                start1 = rawvec.index(['COD1']) +1
                stop1 = rawvec.index(['COD2']) -1
                COD1 = (list(flatten(rawvec[start1: stop1 ])))
                start2 = 1 + rawvec.index(['NONC'])
                NONC = (list(flatten(rawvec[start2: start2 +64 ])))
                itemvect = COD1 + NONC
                #itemvect = [float(i) for i in itemvect]
                assert len(itemvect) == 256, "a full length vector could not be extracted from %s" % f.name
                modvect =  taxlist + itemvect
            f.close()       
            return modvect
            break
            

# for each sequence in the fasta file:
for record in records:
    if not os.path.isdir(path):
        os.mkdir(path) # Generate a temp dir
    os.chdir(path) 
    handle = open("fragment.fasta", "w") # open a fasta file
    SeqIO.write(record, handle, "fasta") # write the sequence to it
    handle.close() # close the file
    
    ## Run metamark
    metamarkparams = ["gmsn.pl", "--clean", "--gm", "--par", mmp,"fragment.fasta"]
    p1 = subprocess.Popen(metamarkparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    metamarkout, metamarkerr= p1.communicate()
    error = re.search("ERROR", metamarkout)
    if error:
        os.chdir(os.path.abspath(args.tmpdir))
        shutil.rmtree(path)
        break
    #Parse results
    featurevect = parsemod(path,args.taxid)
    if featurevect:
        args.output.write("\t".join(featurevect))
        args.output.write("\n")
        os.chdir(os.path.abspath(args.tmpdir))
        shutil.rmtree(path)
    
