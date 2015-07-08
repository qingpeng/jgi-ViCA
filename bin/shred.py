#!/usr/bin/env python

import os
import argparse
import subprocess
import numpy as np
from Bio import SeqIO, SeqRecord

parser = argparse.ArgumentParser(description='A script to shred a genome file into sizes specified by a fixed length or a gamma distbution')
parser.add_argument('--shred', help="Select method to shred the genome contigs", choices =["gamma", "fixed"], default="fixed")
parser.add_argument('--input', help="A fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help= "Name of output file written", type=argparse.FileType('w'), default='-')
parser.add_argument('--samples', help="Total number of shreded contigs to create, or if between 0 and 1, the proportion of the genome to sample", default = 0.5, type=float) 
parser.add_argument('--length', help="The length of the genome subsamples", default = 100, type=int)
parser.add_argument('--shape', help="Shape parameter of gamma distribution, k", default =100, type=float)
parser.add_argument('--scale', help="Scale parameter of gamma distribution, theta", default = 10, type=float)
parser.add_argument('--offset', help= "Offset, or minimum contig length allowed", default = 0, type=int)

args = parser.parse_args()

#Functions
def estimate_samples_gamma(genomelength, samples, shape, scale, offset):
    """Return the number of sequences to output from a gamma distribution"""
    assert samples > 0 , "Sample value must be greater than 0"
    if samples <  1:
        num = genomelength * samples /(shape * scale +offset)
    else: 
        num = samples
    return num

def estimate_samples_fixed(genomelength, samples, length):
    """Return the number of sequences to output from a fixed length"""
    assert samples > 0 , "Sample value must be greater than 0"
    if samples <  1:
        num = genomelength * samples /length
    else: 
        num = samples
    return num

def shred(fasta, shred, samples, shape, scale, offset, length):
    """Take a large FASTA and Return a Multi sequence FASTA with fragments of specified length"""
    #create weighting vectors based on the length of the contigs in the genome
    records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    sampled_frags = []
    slen = [] # the length of contigs in a genome
    ids = [] # the ids of contigs in a genome
    for id,record in records.items():
        ids.append(id)  # add contig ids
        slen.append(len(record.seq)) # add contig lengths to a list
    totlen = sum(slen) # total length of all contigs
    assert totlen > length, "Plese choose a shorter sample length"
    weights = []
    for x in slen:
        weights.append(float(x) / totlen) # Weighting vector for selecting contigs
    #Determine the number of reads to sample
    if shred == "gamma":
        numsamples = estimate_samples_gamma(sum(slen), samples, shape, scale, offset)
    if shred == "fixed":
        numsamples = estimate_samples_fixed(sum(slen), samples, length)
    #select the contigs
    reads_selected = 0
    attempts = 0
    while reads_selected < numsamples:
        attempts += 1
        assert (attempts < 10 * numsamples), "Too many attempts were made subsampling the genome, adjust the sampling parameters"
        #Calculate sample read length based on gamma distribution
        if shred == "gamma":
            length = int(np.random.gamma(shape=shape,scale=scale)) + int(offset)
        print(length)
        id = np.random.choice(a=ids,p=weights)
        selectedseq = records[id]
        selectlen = len(selectedseq)
        try:
            maxstart = selectlen - length
            startpos = np.random.choice(maxstart)
            endpos = startpos+length
            subrecord = selectedseq[startpos:endpos]
            subid = str(id) + str(startpos) + ".." + str(endpos)
            subrecord.id = subid
            sampled_frags.append(subrecord)
            reads_selected += 1
        except:
            pass
    return sampled_frags

samples_frags = shred(fasta=args.input, shred=args.shred, samples=args.samples, \
shape=args.shape, scale=args.scale, offset=args.offset, length=args.length)
SeqIO.write(samples_frags, args.output, "fasta")