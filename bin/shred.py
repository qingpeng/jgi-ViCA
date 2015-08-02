#!/usr/bin/env python

import os
import argparse
import subprocess
import numpy as np
import scipy.stats
from Bio import SeqIO, SeqRecord

parser = argparse.ArgumentParser(description='A script to shred a genome file into sizes specified by a fixed length or a lognormal distbution')
parser.add_argument('--shred', help="Select method to shred the genome contigs", choices =["lognorm", "fixed"], default="fixed")
parser.add_argument('--input', help="A fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help= "Name of output file written", type=argparse.FileType('w'), default='-')
parser.add_argument('--samples', help="Total number of shreded contigs to create, or if between 0 and 1, the proportion of the genome to sample", default = 0.5, type=float) 
parser.add_argument('--length', help="The length of the genome subsamples if fixed is selected", default = 5000, type=int)
parser.add_argument('--shape', help="Shape parameter of lognormal distribution, size", default =1.333, type=float)
parser.add_argument('--loc', help= "Offset, or minimum contig length allowed determind from fitting", default = 3000, type=int)
parser.add_argument('--scale', help="Scale parameter of gamma distribution, theta", default = 1140, type=float)


args = parser.parse_args()

## Functions
def estimate_samples_lognorm(genomelength, samples, shape, loc, scale):
    """Return the number of sequences to output from a gamma distribution"""
    assert samples > 0 , "Sample value must be greater than 0"
    if samples <  1:
    	mean = scipy.stats.lognorm.stats(s = shape, loc = loc, scale = scale, moments='m')
        num = genomelength * samples /mean
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

def shred(fasta, shred, samples, shape, loc, scale, length):
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
    weights = []
    for x in slen:
        weights.append(float(x) / totlen) # Weighting vector for selecting contigs
    #Determine the number of reads to sample
    if shred == "lognorm":
        numsamples = estimate_samples_lognorm(genomelength = totlen, samples = samples,  shape = shape, loc = loc, scale = scale )
    if shred == "fixed":
        numsamples = estimate_samples_fixed(genomelength = totlen, samples = samples, length = length)
    #select the contigs
    reads_selected = 0
    attempts = 0
    while len(sampled_frags) < numsamples:
        attempts += 1
        assert (attempts < 10 * numsamples), "Too many attempts were made subsampling the genome, adjust the sampling parameters"
        #Calculate sample read length based on gamma distribution
        if shred == "lognorm":
            length = int(scipy.stats.lognorm.rvs(s = shape, loc = loc, scale = scale))
        id = np.random.choice(a=ids,p=weights)
        selectedseq = records[id]
        selectlen = len(selectedseq)
        try:
            maxstart = selectlen - length
            startpos = np.random.choice(maxstart)
            endpos = startpos+length
            subrecord = selectedseq[startpos:endpos]
            subid = str(id) + "|" + str(startpos) + ".." + str(endpos)
            subrecord.id = subid
            sampled_frags.append(subrecord)
        except:
            continue
    return sampled_frags


## main

samples_frags = shred(fasta = args.input, shred = args.shred, samples = args.samples, \
shape = args.shape, loc=args.loc, scale = args.scale, length=args.length)
if samples_frags:
	SeqIO.write(samples_frags, args.output, "fasta")