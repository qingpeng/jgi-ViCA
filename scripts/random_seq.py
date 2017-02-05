from Bio import SeqIO
import random

virus_num = float(5000)/963200
nonvirus_num = float(5000)/7745400


input_seq_iterator = SeqIO.parse("virus_segment.fa", "fasta")
virus_seq_iterator = (record for record in input_seq_iterator \
                      if random.random() <=virus_num)


SeqIO.write(virus_seq_iterator, "virus_segment_5k.fa", "fasta")

input_seq_iterator = SeqIO.parse("nonvirus_segment.fa", "fasta")

nonvirus_seq_iterator = (record for record in input_seq_iterator \
                      if random.random() <=nonvirus_num)


SeqIO.write(nonvirus_seq_iterator, "nonvirus_segment_5k.fa", "fasta")

