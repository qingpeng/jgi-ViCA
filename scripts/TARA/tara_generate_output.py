import sys
import random


vfam_obj = open("all_2k_contigs.fa.vfam", 'r')
pfam_obj = open("all_2k_contigs.fa.protein.pfam", 'r')
hit_vfam = set()
hit_pfam = set()

diamond_obj = open("all_2k_contigs.fa.diamond.virus", 'r')

hit_diamond = set()

for line in diamond_obj:
    line = line.rstrip()
    fields = line.split()
    hit_diamond.add(fields[0])
    print fields


for line in vfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        fields = line.split()
        hit_vfam.add(fields[2])

for line in pfam_obj:
    line = line.rstrip()
    if line[0] != '#':
        # print line
        fields = line.split()
        #print fields[4], fields[-1]
        # print fields
        if float(fields[4]) <= 0.00001:

            if ('RNA dependent RNA polymerase' in line or
                'Minor capsid protein' in line or
                    'Major capsid protein' in line):

                hit_pfam.add(fields[2][:-2])

#print hit_pfam
prediction_obj = open("all_2k_prediction.out", 'r')
output_obj = open("all_2k_prediction.out.label", 'w')

for line in prediction_obj:
    line = line.rstrip()

    fields = line.split()
    if fields[0] in hit_vfam:
        line = line + ' ' + '1.0'
    else:
        line = line + ' ' + '0.0'
    if fields[0] in hit_pfam:
        line = line + ' ' + '1.0'
    else:
        line = line + ' ' + '0.0'
    if fields[0] in hit_diamond:
        line = line + ' ' + '1.0'
    else:
        line = line + ' ' + '0.0'

    if fields[0] in hit_vfam or fields[0] in hit_pfam or fields[0] in hit_diamond:
        line = line + ' ' + '1.0'
    else:
        line = line + ' ' + '0.0'


    output_obj.write(line+'\n')
