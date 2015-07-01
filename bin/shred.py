#!/usr/bin/env python

import os
import argparse
import subprocess
import numpy as np
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to shred a genome file into sizes specified by a fixed length or a gamma distbution')
parser.add_argument('--shread', help="Select to shred the genome contigs", choises =["gamma", "fixed","unaltered"], default="fixed")
parser.add_argument('--input', type=argparse.FileType('r'), help="A fasta file", default='-')
parser.add_argument('--output', help= "name of output file written to path", type=argparse.FileType('w'), default='-')
parser.add_argument('--samples', help="Total number of shreded contigs to create, if between 0 and 1 the proportion of the genome to sample", default = 0.5, type=float) 
parser.add_argument('--length', help="The length of the genome subsamples", default = 3000, type=int)
parser.add_argument('--shape', help="Shape parameter of gamma distribution, k", default = 1, type=float)
parser.add_argument('--scale', help="Scale parameter of gamma distribution, theta", default = 1000, type=float)
parser.add_argument('--offset', help= "offset, or minimum contig length allowed", default = 0, type=float)

args = parser.parse_args()

#Functions

def estimate_samples_gamma(genomelength, samples, shape, scale, offset):
    """Return the number of sequences to output from a gamma distribution"""
    assert samples > 0 , "Sample value must be greater than 0"
    if sample <=  1:
        num = genomelength * samples /(shape * scale +offset)
    else: 
        num = samples
    return num

def estimate_samples_fixed(genomelength, samples, length):
    """Return the number of sequences to output from a fixed length"""
    assert samples > 0 , "Sample value must be greater than 0"
    if sample <=  1:
        num = genomelength * samples /length
    else: 
        num = samples
    return num

def shred(fasta, shread, samples, shape, scale, offset, length):
    """Take genome fasta and shred it into pieces"""
    #create weighting vectors based on the length of the contigs in the genome
    records = SeqIO.index(fasta, "fasta")
    sampled_frags = []
    slen = [] # the length of contigs in a genome
    ids = [] # the ids of contigs in a genome
    for record in records:
        ids.append(record.id)  # add contig ids
        slen.append(len(record.seq)) # add contig lengths
    weights = slen/sum(slen) # Weighting vector for selecting contigs
    
    #Determine the number of reads to sample
    if shread == "gamma":
    	numsamples = estimate_samples_gamma(sum(slen), samples, shape, scale, offset)
    if shread == "fixed":
    	numsamples = estimate_samples_fixed(sum(slen), samples, length)
    #select the contigs
    reads_selected = 0
    attempts = 0
    while reads_selected < numsamples:
        attempts += 1
        assert (attempts > 10*samples), "Too many attempts were made subsampling the genome"
        #Calculate sample read length based on gamma distribution
        if len is None:
            length = np.random.gamma(shape=shape,scale=scale,size=1).round + offset
        else:
            length = len
        id = np.random.choice(a=ids,p=slen)
        selectedseq = records[id]
        maxstart = len(selectesseq.seq)-length
        if maxstart < length:
            pass
        startpos = np.random.choice(maxstart)
        endpos = startpos+length
        subrecord = record[startpos:endpos]
        ssr=SeqRecord(subrecord,'fragment_%i' % (i+1),'','')
        sampled_frags.append(ssr)
        reads_selected += 1
    return samples_frags


if args.shread == "unaltered":
    args.output.write(args.input)
elif args.shread == "gamma" or "fixed":
	samples_frags = shred(fasta=args.input, shread=args.shread, samples=args.samples, \
	shape=args.shape, scale=args.scale, offset=args.offset, length=args.length)
	SeqIO.write(samples_frags, args.output, "fasta")