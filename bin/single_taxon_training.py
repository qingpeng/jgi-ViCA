#!/usr/common/usg/languages/python/2.7.4/bin/python

import os
import argparse
import subprocess
import numpy as np

parser = argparse.ArgumentParser(description='A script to run genmark on a single taxon optionally shredding the taxon into sizes spefified by a distbution and creating a file with vector(s) describing the taxon')
parser.add_argument('--path',help="Path containing the sample fasta file")
parser.add_argument('--output', help= "name of output file written to path")
parser.add_argument('--shread', help="Select to shred the genome contigs" action="store_true")
parser.add_argument('--shape', help="Shape parameter of gamma distribution", type=float)
parser.add_argument('--scale', help="scale parameter of gamma distribution", type=float)
parser.add_argument('--samples', help="Total number of shreded contigs to create", type=float) 
parser.add_argument('--ofset', help= "offset, or minimum contig length allowed", type=float)

args = parser.parse_args()

if args.shread:
	lengths = numpy.random.gamma(shape=args.shape,scale=args.scale,size=args.samples)
	lengths.round(0)
	lengths += args.offset

	


