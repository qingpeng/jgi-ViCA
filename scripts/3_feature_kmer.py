#!/usr/bin/env python
import khmer
import argparse
from Bio import SeqIO
from itertools import chain
import itertools
from Bio.Seq import Seq


def iterate_kmer(k):
    """ get the list of tetramers"""
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
        
        
def get_composition(ksize, seq, kmers, norm):
    """ get the composition profile, add one extra count to avoid 0 count"""
    counting_hash = khmer.new_counting_hash(ksize, 200000, 1)
    counting_hash.consume(seq)
    composition = [counting_hash.get(kmer)+1 for kmer in kmers]
    if norm == True:
        total = sum(composition)
        composition_norm = [str(number*1.0/total) for number in composition]
        composition = composition_norm
    return zip(kmers,composition)
    
def generate_line(zip_list):
    """ generate line to printout from list of tuples with feature_name and value"""
    line = ''
    for tuple in zip_list:
        pair = tuple[0]+":"+ str(tuple[1])
        line = line + ' ' + pair
    return line[1:]
    


def main():

    parser = argparse.ArgumentParser(description='A script to generate k-mer coposition frequency')
    parser.add_argument('--input', help="A multi-sequence fasta file", default='-')
    parser.add_argument('--output', help= "Output file, space delimited format", default='-')
    parser.add_argument('--ksize', help="size of kmer, default=4", default = 4)

 
    args = parser.parse_args()

    ## File parsing and variable assignment
    ksize = int(args.ksize)
    records =  SeqIO.parse(args.input, "fasta")
   # print len(records)
 #   print "here"
    kmers = iterate_kmer(ksize)
    file_output_obj = open(args.output, 'w')
# Output format:
# seq_id'\t'seq_length'\t'seq_des'\t'vectors, separated by " "
#


    for record in records:
     #   print record
        if record.description == '':
            des = 'DESCRIPTION'
        else:
            des = record.description
            
        length = len(record.seq)
        kmer_frequency = get_composition(ksize,str(record.seq).upper(), kmers,True)
        
        line = record.id + '\t' + str(length)+'\t'+ des+ '\t'+generate_line(kmer_frequency)
        print line
        file_output_obj.write(line+'\n')
    
    file_output_obj.close()
    
    
if __name__ == '__main__':
    main()