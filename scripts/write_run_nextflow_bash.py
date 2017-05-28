#!/usr/bin/env python

import os

raw_dir = '/global/dna/projectdirs/MEP/tararna/Assemblies/'

dir_list = os.listdir(raw_dir)

group = 1
file_out_obj = open('run_nextflow_'+str(group)+'.sh', 'w')
count = 0
for each_dir in dir_list:
    if count == 12:
        file_out_obj.close()
        group += 1
        file_out_obj = open('run_nextflow_' + str(group) + '.sh', 'w')
        count = 0
    command = 'nextflow run feature_extraction.nf --segment_file  ' \
              ''+raw_dir+each_dir+'/rna/transcripts.fasta  --out  ' \
              '/global/projectb/scratch/qpzhang/TARA/'+each_dir+'.vect'

    file_out_obj.write(command+'\n')
    count += 1
