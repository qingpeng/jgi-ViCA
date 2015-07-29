#!/usr/bin/env python
import argparse
import os
import subprocess
from Bio import SeqIO
from itertools import chain
import shutil
import re
import sys
import tempfile


parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
using emmission data from Metamark')
parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--outfile', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
parser.add_argument('--taxid', help="The taxonomy id")
parser.add_argument('--label', help="Choice of label, normally taxid, but readid for bining applications", choices=['taxid','readid'],default='taxid')
parser.add_argument('--mmp', help="the parameters file for metamark", default = "../gm_parameters/par_11.modified")
args = parser.parse_args()

## File parsing and variable assignment

mmp = os.path.abspath(args.mmp)

# Functions
def flatten(listOfLists):
    """Flatten a list to one level of nesting"""
    return chain.from_iterable(listOfLists)

def parsemod(dir, id):
    """Finds the Genemark model, parses it and returns a feature vector with the taxid as the first element"""
    for file in os.listdir(dir):
        if file.endswith("hmm.mod"):
            with open(file, "r") as f:
                rawvec = []
                idlist = [id]
                for line in f.readlines():
                    rawvec.append(line.split())
                start1 = rawvec.index(['COD1']) +1
                stop1 = rawvec.index(['COD2']) -1
                COD1 = (list(flatten(rawvec[start1: stop1 ])))
                start2 = 1 + rawvec.index(['NONC'])
                NONC = (list(flatten(rawvec[start2: start2 +64 ])))
                itemvect = COD1 + NONC
                if len(itemvect) == 256:
                    modvect =  idlist + itemvect
                    return modvect
            f.close()
            break
            

# for each sequence in the fasta file:
records = SeqIO.parse(args.input, "fasta")
cnt_success = 0
cnt_vectfailure = 0
cnt_mmfailure = 0
len_records =0
for record in records:
    len_records += 1
    tmpdir = tempfile.mkdtemp(dir="/scratch")
    os.chdir(tmpdir) 
    handle = open("fragment.fasta", "w") # open a fasta file
    SeqIO.write(record, handle, "fasta") # write the sequence to it
    handle.close() # close the file
    ## Run metamark
    metamarkparams = ["gmsn.pl", "--clean", "--gm", "--par", mmp,"fragment.fasta"]
    p1 = subprocess.Popen(metamarkparams, stdout=subprocess.PIPE)
    metamarkout, metamarkerr= p1.communicate()
    if args.label == 'taxid':
    	featurevect = parsemod(tmpdir,args.taxid)
    elif args.label == 'readid':
    	featurevect = parsemod(tmpdir,record.id)
    else:
    	raise InputError("the label parameter must be either 'taxid' or 'readid'")
    if p1.returncode == 0:
        featurevect = parsemod(tmpdir,args.taxid)
        if featurevect:
            args.outfile.write("\t".join(featurevect))
            args.outfile.write("\n")
            shutil.rmtree(tmpdir)
            cnt_success += 1
        else:
            shutil.rmtree(tmpdir)
            cnt_vectfailure  += 1
    else:
        shutil.rmtree(tmpdir)
        cnt_mmfailure += 1

if cnt_success == 0:
	args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, metamark errors: %s, vector errors: %s \n" \
	% (args.taxid, len_records, cnt_success, cnt_mmfailure, cnt_vectfailure))
args.outfile.close()

