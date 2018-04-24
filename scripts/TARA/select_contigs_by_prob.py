import sys
import random
from Bio import SeqIO
import sys
import os

file_in = "all_4k_prediction.out"

result = open(file_in, 'r')

contig_list = set()
for line in result:
    line = line.rstrip()
    fields = line.split()
    if float(fields[2]) >= 0.99995837748:
        contig_list.add(fields[0])


input_seq_iterator = SeqIO.parse(open("all_4k_contigs.fa", "rU"), "fasta")
long_list = []
for record in input_seq_iterator:
 #   print record.name
 #   print record.id
 #   exit()
    if record.name in contig_list:
        long_list.append(record)



output_handle = open("all_4k_contigs_top1k.fa", "w")
SeqIO.write(long_list, output_handle, "fasta")
output_handle.close()