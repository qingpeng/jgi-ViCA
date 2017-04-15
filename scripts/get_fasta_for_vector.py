#!/usr/bin/env python

from Bio import SeqIO


file_vector_list = ['non_virus_testing_5k.vect', 'non_virus_training_25k.vect',
                    'virus_testing_5k.vect', 'virus_training_25k.vect']

file_fasta = 'all_segment.fa'

big_list = []
seg_name_set = set()

for file_vector in file_vector_list:
    file_vector_obj = open(file_vector, 'r')
    seg_name_list = []
    for line in file_vector_obj:
        line = line.rstrip()
        fields = line.split()
        seg_name = fields[0]
        seg_name_list.append(seg_name)
        seg_name_set.add(seg_name)
    big_list.append(seg_name_list)
    print file_vector+'done!\n'


input_seq_iterator = SeqIO.parse(file_fasta, "fasta")

fasta_dict = {}
for record in input_seq_iterator:
    if record.id in seg_name_set:
        fasta_dict[record.id] = record
print 'all.fasta done!\n'


for i in range(4):
    file_out = file_vector_list[i]+'.fasta'
    records = []
    for seg_name in big_list[i]:
        records.append(fasta_dict[seg_name])

    SeqIO.write(records, file_out, "fasta")
    print file_out+' done!\n'
