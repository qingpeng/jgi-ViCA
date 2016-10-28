#! /usr/bin/env python
"""
Generate 4-mer frequency profile for each segment.

% feature_extraction_kmer.py <segment file> <profile output>

Use '-h' for parameter help.
"""

import sys
import os
# require khmer 1.4.1
import khmer
import argparse
import itertools
from Bio.Seq import Seq
from Bio import SeqIO

def iterate_kmer(k):
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
    counting_hash = khmer.new_counting_hash(4, 2000, 1)
    counting_hash.consume(seq)
    composition = [counting_hash.get(kmer) for kmer in kmers]
    if norm == True:
        total = sum(composition)
        composition_norm = [str(number*1.0/total) for number in composition]
        composition = composition_norm
    return composition
    

def get_composition2(seq, kmers, norm):
    composition_hash = {}
    for i in range(len(seq)-3):
        try:
            composition_hash[seq[i:i+4]] += 1
        except:
            composition_hash[seq[i:i+4]] = 1
#    print composition_hash
    composition = []
    
    str(Seq(seq[i:i+4]).reverse_complement())
    for kmer in kmers:
        try:
            both_count = composition_hash[kmer]
        except:
            both_count = 0
        if kmer != str(Seq(kmer).reverse_complement()):
            try:
                both_count = both_count 
                + composition_hash[str(Seq(kmer).reverse_complement())]
            except:
                both_count = both_count
        composition.append(both_count)

    if norm == True:
        total = sum(composition)
        composition_norm = [number*1.0/total for number in composition]
        composition = composition_norm
    return composition
    
def main():


    parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
    using k-mer frequency profile')
    parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
    parser.add_argument('--outfile', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
    parser.add_argument('--taxid', help="The taxonomy id")
    parser.add_argument('--label', help="Choice of label, normally taxid, but readid for bining applications", choices=['taxid','readid'],default='taxid')
    parser.add_argument('--minlen', help="minimum length to attempt to classify", default = 3000)
    args = parser.parse_args()
    
    records = SeqIO.parse(args.input, "fasta")
    kmers = iterate_kmer(4)
    shortreads = 0
    cnt_success = 0
    len_records =0
    for record in records:
        len_records = len_records + 1
        if len(record) < args.minlen:
            shortreads += 1
            continue
        
        norm = True
        if args.label == 'taxid':
            vect = [args.taxid] + [record.description] + get_composition(str(record.seq).upper(), kmers,norm)
            args.outfile.write("\t".join(vect))
            args.outfile.write("\n")
        elif args.label == 'readid':
            vect = [record.id] + [record.description] + get_composition(str(record.seq).upper(), kmers,norm)
            args.outfile.write("\t".join(vect))
            args.outfile.write("\n")
        cnt_success += 1
    if cnt_success == 0:
        #args.outfile.write("#Taxon id: 123, Number of Contigs: 3\n") #kmer_v3
        args.outfile.write("#Taxon id: 123\n") #kmer_v4
       # args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, reads below metamark min: %s \n" \
        #% (args.taxid, len_records, cnt_success, shortreads))
   # print "#Taxon id:",args.taxid,"Sucessful contigs:",len_records,"Failed contigs:",shortreads,"\n"
    args.outfile.close()
    
if __name__ == '__main__':
    main()

