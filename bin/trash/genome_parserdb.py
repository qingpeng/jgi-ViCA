#!/usr/common/usg/languages/python/2.7.4/bin/python
#genome_parserdb.py
#Adam Rivers 2015/03/19 arrivers@lbl.gov
import argparse
import gzip
import os.path
import re
import sqlite3
from Bio import SeqIO
import sys

parser = argparse.ArgumentParser(description='A script to split NCBI RefSeq genome archive files into fasta files containing a single genome per file')
parser.add_argument('-i','--input',help='A gzip compressed fasta input file or an uncompressed stdin if -')
parser.add_argument('-o','--out', help='A user supplied output directory for genome files', required =True)
parser.add_argument('-t','--taxdb', help='The location of the NCBI taxonomy database', required =True)
args = parser.parse_args()

#parse stdin vs file
if args.input == '-':
	refseqfile = sys.stdin
else:
	refseqfile = gzip.open(args.input, 'rb')

# output path
abspath = os.path.abspath(args.out)

#SQLite Database Connection
conn = sqlite3.connect(os.path.abspath(args.taxdb))
c = conn.cursor()

# Create dictionary to store taxa that have been seen
taxadict = {}

#Fasta Record look up the taxonomy and place into a genome file
for record in SeqIO.parse(refseqfile, "fasta"):
	gi = re.sub(r'\|.*','',re.sub(r'gi\|','',record.id)) # get gi
	git = (gi,) #make gi into tuple
	c.execute('SELECT taxid from gi_taxid_nucl where gi=?',git) # query DB for taxid
	taxid = c.fetchone() # retrieve answer
	if taxid is not None: # skip any sequences with no taxonomy information
		if taxid in taxadict: # write sequences to existing genome files if present
			fname = abspath + "/"+ str(taxid[0]) + ".fasta"
			with open(fname, 'a') as genomefile:
				SeqIO.write(record, genomefile, "fasta")
		else: # write sequences to new files if no other sequences are yet present from the genome
			fname = abspath + "/"+ str(taxid[0]) + ".fasta"
			with open(fname, 'a') as genomefile:
				SeqIO.write(record, genomefile, "fasta")
				taxadict[taxid[0]] = None # add taxonomy ID to the dictionary
refseqfile.close()