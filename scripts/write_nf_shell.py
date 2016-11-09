#!/usr/bin/env python

for i in range(20):
    if i <10:
        file_name = "segment_all.000"+str(i)+".fa"
    else:
        file_name = "segment_all.00"+str(i)+".fa"

    cmd = "nextflow run ~/Bitbucket/jgi-genelearn/scripts/get_vector_lite.nf -qs 50 --out=./out"+str(i)+".vect --segment_file=./Split/"+file_name + ' &'
    print cmd

