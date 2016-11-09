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



# Functions
def flatten(listOfLists):
    """Flatten a list to one level of nesting"""
    return chain.from_iterable(listOfLists)

def parsemod(dir):
    """Finds the Genemark model, parses it and returns a feature vector with the taxid as the first element"""
    for file in os.listdir(dir):
        if file.endswith("hmm.mod"):
            file = os.path.join(dir,file)
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
                if len(itemvect) == 256:
                    return itemvect
            f.close()
        break



def main():

    parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
    using emmission data from GeneMarkS')
    parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
    parser.add_argument('--outfile', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
    parser.add_argument('--taxid', help="The taxonomy id")
    parser.add_argument('--label', help="Choice of label, normally taxid, but readid for binning applications", choices=['taxid','readid'],default='taxid')
    parser.add_argument('--mmp', help="the parameters file for genemark", default = "../gm_parameters/par_11.modified")
    parser.add_argument('--tmp', help="root directory to write temp files in", default = "/scratch")
    parser.add_argument('--minlen', help="minimum length to attempt to classify", default = 3000)
    parser.add_argument('--prog', help="metamerk program to run", choices=['genemarks','metagenemark'],default='genemarkS')
    args = parser.parse_args()

    ## File parsing and variable assignment

    mmp = os.path.abspath(args.mmp)
    tmp = os.path.abspath(args.tmp)

    # for each sequence in the fasta file:
    records = SeqIO.parse(args.input, "fasta")
    cnt_success = 0
    cnt_vectfailure = 0
    cnt_mmfailure = 0
    len_records =0
    shortreads = 0

    for record in records:
        #go on if reads are too short
        if len(record) < args.minlen:
            shortreads += 1
            continue
        len_records += 1

        tmpdir = tempfile.mkdtemp(dir=tmp)
        os.chdir(tmpdir) 
        handle = open("fragment.fasta", "w") # open a fasta file
        SeqIO.write(record, handle, "fasta") # write the sequence to it
        handle.close() # close the file
        
        ## Run genemarks
        genemarkparams = ["gmsn.pl", "--clean", "--gm", "--par", mmp,"fragment.fasta"]
        p1 = subprocess.Popen(genemarkparams, stdout=subprocess.PIPE)
        genemarkout, genemarkerr= p1.communicate()
        if p1.returncode == 0:
            featurevect = parsemod(tmpdir)
            if featurevect:
                if args.label == 'taxid':
                    vect = [args.taxid] + [record.description] + featurevect
                elif args.label == 'readid':
                    vect = [record.id] + [record.description] + featurevect
                else:
                    raise InputError("the label parameter must be either 'taxid' or 'readid'")
                args.outfile.write("\t".join(vect))
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
        args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, genemark errors: %s, vector errors: %s, reads below Genemark min: %s \n" \
        % (args.taxid, len_records, cnt_success, cnt_mmfailure, cnt_vectfailure,shortreads))
args.outfile.close()
    
if __name__ == '__main__':
    main()
