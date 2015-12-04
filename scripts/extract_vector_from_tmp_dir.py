#!/usr/bin/env python
import argparse
import os



parser = argparse.ArgumentParser(description='A script to extract genomic features from TMP directory')
parser.add_argument('-d', '--directory', help ='TMP directory',required=True) 
parser.add_argument('-o', '--outfile', help ='output file', required=True) 

args = parser.parse_args()



def parse_vect(dir, output):
    """parse vectors from tmp directory"""
    fh_output = open(output, 'w')
    
    for dir2 in os.listdir(dir):
        print dir2
        if os.path.isdir(os.path.join(dir,dir2)):
            for file in os.listdir(os.path.join(dir,dir2)):
                file_name = os.path.join(dir,dir2,file)
                print file_name
                f = open(file_name, 'r')
                for line in f:
                    if line[0] != '#':
                        fh_output.write(line)
    fh_output.close()
    

parse_vect(args.directory, args.outfile)


