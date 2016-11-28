#!/usr/bin/env python

# using full genome, no shredding
# genemarkS + tetramer profile




import argparse
import os
import subprocess
from Bio import SeqIO
from itertools import chain
import shutil
import re
import sys
import tempfile
# require khmer 1.4.1
import khmer
from Bio.Seq import Seq
import itertools



def iterate_kmer(k):
    """ get the list of tetramers"""
    bases = ['A','C','T','G']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
#    print kmers
    core_kmer = []
    for kmer in kmers:
        if not str(Seq(kmer).reverse_complement()) in core_kmer:
            core_kmer.append(kmer)
#    print core_kmer
#    print len(core_kmer)
    return core_kmer
        
        
def get_composition(seq, kmers, norm):
    """ get the composition profile, add one extra count to avoid 0 count"""
    counting_hash = khmer.new_counting_hash(4, 2000, 1)
    counting_hash.consume(seq)
    composition = [counting_hash.get(kmer)+1 for kmer in kmers]
    if norm == True:
        total = sum(composition)
        composition_norm = [str(number*1.0/total) for number in composition]
        composition = composition_norm
    return composition
    
    
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
    using emmission data from Metamark')
    parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
    parser.add_argument('--outfile', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
    parser.add_argument('--taxid', help="The taxonomy id")
    parser.add_argument('--label', help="Choice of label, normally taxid, but readid for bining applications", choices=['taxid','readid'],default='taxid')
    parser.add_argument('--mmp', help="the parameters file for metamark", default = "../gm_parameters/par_11.modified")
    parser.add_argument('--tmp', help="root directory to write temp files in", default = "/scratch")
    parser.add_argument('--minlen', help="minimum length to attempt to classify", default = 3000)
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
    kmers = iterate_kmer(4)
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
        ## Run metamark
        metamarkparams = ["gmsn.pl", "--clean", "--gm", "--par", mmp,"fragment.fasta"]
        p1 = subprocess.Popen(metamarkparams, stdout=subprocess.PIPE)
        metamarkout, metamarkerr= p1.communicate()
        if p1.returncode == 0:
            featurevect = parsemod(tmpdir)
            if featurevect:
                if args.label == 'taxid':
                    vect = [args.taxid] + [record.description] + featurevect + \
                    get_composition(str(record.seq).upper(), kmers,True)
                elif args.label == 'readid':
                    vect = [readid] + [record.description] + featurevect + \
                    get_composition(str(record.seq).upper(), kmers,True)
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
        args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, metamark errors: %s, vector errors: %s, reads below metamark min: %s \n" \
        % (args.taxid, len_records, cnt_success, cnt_mmfailure, cnt_vectfailure,shortreads))
    args.outfile.close()
    
if __name__ == '__main__':
    main()