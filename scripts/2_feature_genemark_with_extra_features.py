#!/usr/bin/env python


# feature:
# 3.  get lowest e-value for pfam and vfam, if there are multiple hits,
# move such functions from convert to libsvm to here (09/23 updated)
# 
# hmm profile needs to be compressed before running HMMER
# hmmpress vFam_IMG.hmm
# require module load hmmer/3.1b2
#

import khmer
import argparse
from Bio import SeqIO
from itertools import chain
import itertools
from Bio.Seq import Seq
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

    itemvect = []
    name_list = []
    value_list = []
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

            f.close()

        elif file_name.endswith("fasta.lst"):
            lst_file = os.path.join(dir_name, file_name)
            lst_file_obj = open(lst_file, 'r')
            switch = 0
            gene_lines = []
            segment_length = 0
            for line in lst_file_obj:
                if "FASTA definition line" in line:
                    line = line.rstrip()
                    start_end = line.split()[3].split('|')[-1].split('..')
                    segment_length = int(start_end[1]) - int(start_end[0])

                if "Predicted genes" in line:
                    switch = 1
                elif "Predicted proteins:" in line:
                    switch = 0
                else:
                    if switch == 1:
                        if (('Gene' not in line) and ('#' not in line) and
                                line != '\n'):
                            gene_lines.append(line)
            lst_file_obj.close()

            left_end = []
            right_end = []
            gene_length = []
            # print gene_lines
            for gene_line in gene_lines:
                gene_line = gene_line.rstrip()
                fields = gene_line.split()
                # print fields
                left_end.append(int(fields[2].strip('<')))
                right_end.append(int(fields[3].strip('>')))
                gene_length.append(int(fields[4]))

            fraction_coding = sum(gene_length)*1.0/segment_length
            ave_gene_length = sum(gene_length)*1.0/len(gene_length)

            if len(gene_length) == 1:
                ave_intergenic_length = 0
            else:
                ave_intergenic_length = (1.0*((right_end[-1]-left_end[0]+1)
                                         - sum(gene_length)) /
                                         (len(gene_length)-1))
            name_list = ['fraction_coding', 'ave_gene_length',
                         'ave_intergenic_length']
            value_list = [fraction_coding, ave_gene_length,
                          ave_intergenic_length]
            # print value_list, gene_length, segment_length

    id_list = range(256) + name_list
    num_list = itemvect + value_list

    if len(num_list) == 259:
        return zip(id_list, num_list)





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
        metamarkout, metamarkerr= p1.communicate()
        
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
                metamarkout, metamarkerr= p2.communicate()

                if p2.returncode == 0:  # if probuild can generate hmm model

                    vect_genemark = parsemod(tmpdir)
                    line_genemark = (line_genemark +
                                     generate_line(vect_genemark))


            
        output_genemark_obj.write(line_genemark+'\n')

        shutil.rmtree(tmpdir)
    
    output_genemark_obj.close()

    
if __name__ == '__main__':
    main()
