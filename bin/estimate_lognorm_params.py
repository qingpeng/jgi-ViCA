#!/usr/bin/env python

import argparse
import os
import subprocess
from Bio import SeqIO
import sys
import numpy as np
import scipy
import scipy.stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser( \
	description='Estimate-Log normal distribution parameters describing the length of assemblies by sampling')
parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
parser.add_argument('--output', help="output file with the 3 parameters" , type=argparse.FileType('w'), default='-')
parser.add_argument('--fitplot', help ="a file with more information on the fit", default="lognormfit.pdf") 
parser.add_argument('--plotformat', help ="File format to plot the fit graph", choices = ["png", "pdf", "ps", "eps", "svg"], default="pdf") 
args = parser.parse_args()


cnt = []
sequences = SeqIO.parse(args.input, 'fasta')
for sequence in sequences:
	slen = len(sequence)
	cnt.append(slen)

# convert length data to a numpy array
npcnt = np.array(cnt)
loc = min(cnt)
#fit gamma parameters with scipi 
fit_shape, fit_loc , fit_scale=scipy.stats.lognorm.fit(npcnt, loc = loc)
#print(scipy.stats.lognorm.fit(data = npcnt, loc = loc))
#output the parameters
args.output.write("shape"+"\t"+"loc"+"\t"+"fit_scale"+"\n")
args.output.write(str(fit_shape)+"\t"+str(fit_loc)+"\t"+str(fit_scale)+"\n")

# Plot the results for manual inspection

fig, ax = plt.subplots(1, 1)
# Create curve points for fitted line
x = np.linspace(scipy.stats.lognorm.ppf(0.01, s = fit_shape, loc=loc,scale=fit_scale),\
	scipy.stats.lognorm.ppf(0.99, s = fit_shape, loc = loc, scale = fit_scale), 100)
# plot line
ax.plot(x, scipy.stats.lognorm.pdf(x, s = fit_shape, loc = loc, scale = fit_scale),'r-', lw=2, \
	alpha=0.6, label='lognorm pdf')
# plot histogram of data
ax.hist(npcnt, normed=True, alpha=0.2, bins=200, range=[0, 50000])
# Save figure to file
plt.savefig(args.fitplot, transparent=True, format= args.plotformat)

