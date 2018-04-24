#!/usr/bin/python
from Bio import SeqIO
import sys
import os

#usage: python long.seq.py in.fasta out.fasta 200

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), "fasta")
long_list = []
for record in input_seq_iterator:
    if len(record.seq) >= int(sys.argv[3]):
        long_list.append(record)
    else:
        break


output_handle = open(sys.argv[2], "w")
SeqIO.write(long_list, output_handle, "fasta")
output_handle.close()
