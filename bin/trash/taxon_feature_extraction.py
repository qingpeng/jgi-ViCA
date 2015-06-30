#!/usr/bin/env python

import os
import argparse
import subprocess
import numpy as np
from Bio import SeqIO

parser = argparse.ArgumentParser(description='A script to run genmark on a single taxon optionally shredding the taxon into sizes spefified by a distbution and creating a file with vector(s) describing the taxon')
parser.add_argument('--genome',help="Path containing the sample fasta file", required=True)
parser.add_argument('--output', help= "name of output file written to path", required=True)
parser.add_argument('--length', help="The length of the genome subsamples", default=3000, type=int)
parser.add_argument('--samples', help="Total number of shreded contigs to create, if between 0 and 1 the proportion of the genome to sample", default=0.5, type=float) 
parser.add_argument('--shread', help="Select to shred the genome contigs", action="store_true")
#Needed if shred is selected
parser.add_argument('--shape', help="Shape parameter of gamma distribution", type=float)
parser.add_argument('--scale', help="Scale parameter of gamma distribution", type=float)
parser.add_argument('--offset', help= "offset, or minimum contig length allowed", type=float)

args = parser.parse_args()

#Functions
def shred(fasta, samples, shape, scale, offset):
#create weighting vectors based on the length of the contigs in the genome
    records = SeqIO.index(fasta, "fasta")
    sampled_frags = []
    slen = [] # the length of contigs in a genome
    ids = [] # the length of contigs in a genome
    for record in records:
        ids.append(record.id)  # add contig ids
        slen.append(len(record.seq)) # add contig lengths
    weights = slen/sum(slen) # Weighting vector for selecting contigs
    
    #select the contigs
    reads_selected = 0
    attempts = 0
    while reads_selected < samples:
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

