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
parser.add_argument('--tmp', help="root directory to write temp files in", default = "/scratch")
args = parser.parse_args()

## File parsing and variable assignment

mmp = os.path.abspath(args.mmp)

## Functions
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

def parse_fasta(record):
	"""Parses a fasta record id and description and returns a python list with metdata about that record"""
	# metadata = [id, start,stop,taxid, organelle,code, taxonomic_lineage]
	metadata = []
	if re.search('\|pos\|.',record.id):
		id = re.sub('\|pos\|.','',record.id)
		pos = re.sub('.\|pos\|','',record.id)
		start = re.sub('\.\..','',pos)
		end = re.sub('.\.\.','',pos)
		metadata.append(id)
		metadata.append(start)
		metadata.append(stop)
	else:
		id = record.id
		metadata.append(id)
		metadata.append('')
		metadata.append('')
	
	descriptionlist = record.description.split(',')
	for item in descriptionlist:
		if re.search("taxid=", item):
			taxid = re.sub("taxid=",'',item)
			metadata.append(taxid)
		else:
			metadata.append('')
		if re.search("organelle=", item):
			organelle = re.sub("organelle=",'',item)
			metadata.append(organelle)
		else:
			metadata.append('')
		if re.search("plasmid=", item):
			plasmid = re.sub("plasmid=",'',item)
			metadata.append(plasmid)
		else:
			metadata.append('')
		if re.search("code=", item):
			code = re.sub("code=",'',item)
			metadata.append(code)
		else:
			metadata.append('')	
		if re.search("taxonomy=", item):
			taxonomy = re.sub("taxonomy=",'',item)
			metadata.append(taxonomy)
		else:
			metadata.append('')
	return metadata

# read minimum length in metamark config file
with open(mmp, 'r') as metamarkparams:
	for line in metamarkparams:
		if re.search("--MIN_CONTIG_SIZE",line):
			mcs = re.sub('\n','',re.sub("--MIN_CONTIG_SIZE",'',line)).strip()
metamarkparams.close()

# for each sequence in the fasta file:
records = SeqIO.parse(args.input, "fasta")
cnt_success = 0
cnt_vectfailure = 0
cnt_mmfailure = 0
len_records =0
shortreads = 0
for record in records:
    readid = record.id
    if mcs:
    	if int(mcs) > len(record):
    		shortreads += 1
    		continue
    len_records += 1
    tmpdir = tempfile.mkdtemp(dir=args.tmp)
    os.chdir(tmpdir) 
    handle = open("fragment.fasta", "w") # open a fasta file
    SeqIO.write(record, handle, "fasta") # write the sequence to it
    handle.close() # close the file
    ## Run metamark
    metamarkparams = ["gmsn.pl", "--clean", "--gm", "--par", mmp,"fragment.fasta"]
    p1 = subprocess.Popen(metamarkparams, stdout=subprocess.PIPE)
    metamarkout, metamarkerr= p1.communicate()
    if p1.returncode == 0:
        if args.label == 'taxid':
        	featurevect = parsemod(tmpdir,args.taxid)
    	elif args.label == 'readid':
        	featurevect = parsemod(tmpdir,readid)
    	else:
        	raise InputError("the label parameter must be either 'taxid' or 'readid'")
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
    args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, metamark errors: %s, vector errors: %s, reads below metamark min: %s \n" \
    % (args.taxid, len_records, cnt_success, cnt_mmfailure, cnt_vectfailure,shortreads))
args.outfile.close()

