#!/usr/bin/env python
import argparse
import sys
import operator

parser = argparse.ArgumentParser(description='A script to extract possible virus sequence according to predicted score from classifier' )
parser.add_argument('--seq', help="sequence file in fasta",required = True)
parser.add_argument('--output', help="output file with extracted possible virus sequence",required = True)
parser.add_argument('--number', help ="number/percentage of sequences to extract, if <1, percentage, >=1, number",  required = True) 
parser.add_argument('--score', help ="file with predicted score from Spark", required = True) 
parser.add_argument('--vector', help ="file with vector for each sequence", required = True) 

args = parser.parse_args()

fh_seq = open(args.seq, 'r')
fh_output = open(args.output, 'w')
fh_score = open(args.score,'r')
fh_vector = open(args.vector, 'r')
number = float(args.number)

# get list of all the scores
scores = []
for line in fh_score:
	line = line.rstrip()
	fields = line[1:-1].split(',')
	scores.append(float(fields[0]))

print "score done!"

# get dictionary with seq_name as key, and score as value
seq_score = {} 
k = 0
for line in fh_vector:
	line = line.rstrip()
	fields = line.split()
	seq_score[fields[1]] = scores[k]
	k = k + 1
print "dictionary done!"

# sort seqs by score
sorted_seq_score = sorted(seq_score.items(), key=operator.itemgetter(1), reverse=True)
print "sorting done!"

if number < 1:
	num = int(k * number)
else:
	num = int(number)
# seqs to pick
seq_id = sorted_seq_score[0:num]	
print seq_id
from Bio import SeqIO

fasta_seq = SeqIO.parse(fh_seq, 'fasta')
seqs = {}
for fasta in fasta_seq:
	seqs[fasta.id] = fasta
seq_record_list = []
for id_tuple in sorted_seq_score[0:num]:
	seq_id = id_tuple[0]
	seq_score = id_tuple[1]
	seq_record = seqs[seq_id]
	seq_record.description = str(seq_score)+',length:'+str(len(seq_record))
	seq_record_list.append(seq_record)
print "to write!"

SeqIO.write(seq_record_list, fh_output, "fasta")


