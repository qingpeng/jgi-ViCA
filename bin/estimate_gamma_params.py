#!/usr/bin/env python

import argparse
import os
import subprocess
from Bio import SeqIO
import sys
import numpy as np
import scipy
import scipy.stats
from matplotlib import pyplot as plt


parser = argparse.ArgumentParser(description='Estimate gamma parameters describing the length of assemblies by sampling')
parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help="output file with the 4 parameters" , type=argparse.FileType('w'), default='-')
parser.add_argument('--fitdoc', help ="a file with more information on the fit") 
args = parser.parse_args()


cnt = []
sequences = SeqIO.parse(args.input, 'fasta')
for sequence in sequences:
	slen = len(sequence)
cnt.append(slen)

#TODO fit with scipi script
#find min value 
#subtract min from data and write as offset
#fit with scipi fit gamma
