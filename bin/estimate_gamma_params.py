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


parser = argparse.ArgumentParser( \
	description='Estimate gamma parameters describing the length of assemblies by sampling')
parser.add_argument('--input', help="A multi-sequence fasta file",\
	type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help="output file with the 3 parameters" , \
	type=argparse.FileType('w'), default='-')
parser.add_argument('--fitplot', help ="a file with more information on the fit", \
	defualt="gammafit.pdf") 
parser.add_argument('--plotformat', help ="File format to plot the fit graph", \
	choices = ["png", "pdf", "ps", "eps" "svg"],defualt="pdf") 
args = parser.parse_args()


cnt = []
sequences = SeqIO.parse(args.input, 'fasta')
for sequence in sequences:
	slen = len(sequence)
cnt.append(slen)

# convert length data to a numpy array
npcnt = np.array(cnt)

#fit gamma parameters with scipi 
fit_alpha, fit_loc, fit_scale=scipy.stats.gamma.fit(data)

#output the parameters
args.output.write("alpha","\t","loc","\t","fit_beta","\n")
args.output.write(fit_alpha,"\t",fit_loc,"\t",fit_beta,"\n")

# Plot the results for manual inspection
fig, ax = plt.subplots(1, 1)
# Create curve points for fitted line
x = np.linspace(scipy.stats.gamma.ppf(0.01, fit_alpha, loc=fit_loc,scale=fit_scale),\
	scipy.stats.gamma.ppf(0.99, fit_alpha, loc=fit_loc,scale=fit_scale), 100)
# plot line
ax.plot(x, scipy.stats.gamma.pdf(x, fit_alpha, loc=fit_loc,scale=fit_scale),'r-', lw=2, \
	alpha=0.6, label='gamma pdf')
# plot histogram of data
ax.hist(npcnt, normed=True, alpha=0.2)
# Save figure to file
plt.savefig(args.fitplot, transparent=True, format= args.plotformat)
