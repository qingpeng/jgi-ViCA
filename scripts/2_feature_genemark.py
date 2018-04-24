#!/usr/bin/env python


import argparse
from Bio import SeqIO
from itertools import chain
import subprocess
import tempfile
import os
import shutil


def flatten(list_of_lists):
    """Flatten a list to one level of nesting"""
    return chain.from_iterable(list_of_lists)
    
    
def parsemod(dir_name):
    """Finds the Genemark model, parses it and returns a feature vector
    with the taxid as the first element"""

    for file_name in os.listdir(dir_name):
        if file_name.endswith("hmm.mod"):
            model_file = os.path.join(dir_name, file_name)
            with open(model_file, "r") as f:
                rawvec = []
                for line in f.readlines():
                    rawvec.append(line.split())
                start1 = rawvec.index(['COD1']) + 1
                stop1 = rawvec.index(['COD2']) - 1
                cod1 = (list(flatten(rawvec[start1: stop1])))
                start2 = 1 + rawvec.index(['NONC'])
                nonc = (list(flatten(rawvec[start2: start2 + 64])))
                itemvect = cod1 + nonc
                if len(itemvect) == 256:
                    return zip(range(256), itemvect)
            f.close()
            break


def generate_line(zip_list):
    """ generate line to printout from list of tuples
    with feature_name and value"""
    line = ''
    for tuples in zip_list:
        pair = str(tuples[0])+":" + str(tuples[1])
        line = line + ' ' + pair
    return line[1:]


def main():

    parser = argparse.ArgumentParser(
        description='A script to generate feature vector of genemark')
    parser.add_argument('--input',
                        help="A multi-sequence fasta file", required=True)
    parser.add_argument('--output_prefix',
                        help="prefix of output file, space delimited format, "
                        "*.genemark", required=True)
    parser.add_argument('--tmp',
                        help="tmp directory", default="./")
    parser.add_argument('--genemark_path', help="genemark path, ./gmsuite",
                        required=True)

# Output format:
# seq_id'\t'seq_length'\t'seq_des'\t'vectors, separated by " "
#

    args = parser.parse_args()

    # File parsing and variable assignment
    genemark_path = os.path.abspath(args.genemark_path)
    genelearn_path = os.path.dirname(os.path.realpath(__file__))
    tmp = os.path.abspath(args.tmp)

    records = SeqIO.parse(args.input, "fasta")
    
    output_prefix = args.output_prefix

    output_genemark_obj = open(output_prefix+'.genemark', 'w')

    for record in records:
        tmpdir = tempfile.mkdtemp(dir=tmp)
        os.chdir(tmpdir) 
        handle = open("fragment.fasta", "w")  # open a fasta file
        SeqIO.write(record, handle, "fasta")  # write the sequence to it
        handle.close()  # close the file
        
        length = len(record.seq)

        if record.description == '':
            des = 'DESCRIPTION'
        else:
            des = record.description
            
        line_genemark = record.id+'\t'+str(length)+'\t'+des+'\t'

        meta_genemark_params = [genemark_path+"/gmhmmp", "-m", genemark_path +
                                '/MetaGeneMark_v1.mod', "fragment.fasta", "-a",
                                "-A", "fragment.fasta.aa"]
        p1 = subprocess.Popen(meta_genemark_params, stdout=subprocess.PIPE)
        p1.communicate()
        
        if p1.returncode == 0:
            # if MetaGeneMark can generate the gene prediction
                probuild_params = [genemark_path+"/probuild",
                                   "--par", genelearn_path +
                                   '/gm_parameters/par_11.modified',
                                   "--ORDM", "2", "--order_non", "2",
                                   "--revcomp_non", "1",
                                   "--seq", "fragment.fasta",
                                   "--geneset", "fragment.fasta.lst",
                                   "--mkmod", "hmm.mod"]
                p2 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                p2.communicate()

                if p2.returncode == 0:  # if probuild can generate hmm model

                    vect_genemark = parsemod(tmpdir)
                    line_genemark = (line_genemark +
                                     generate_line(vect_genemark))

        output_genemark_obj.write(line_genemark+'\n')

        shutil.rmtree(tmpdir)
    
    output_genemark_obj.close()

    
if __name__ == '__main__':
    main()
