#!/usr/bin/env python
import argparse
import os
import subprocess
from Bio import SeqIO
from itertools import chain
import shutil
import re
import sys
import tempfile
# require khmer 1.4.1
import khmer
from Bio.Seq import Seq
import itertools



def iterate_kmer(k):
    """ get the list of tetramers"""
    bases = ['A','C','T','G']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
#    print kmers
    core_kmer = []
    for kmer in kmers:
        if not str(Seq(kmer).reverse_complement()) in core_kmer:
            core_kmer.append(kmer)
#    print core_kmer
#    print len(core_kmer)
    return core_kmer
        
        
def get_composition(seq, kmers, norm):
    """ get the composition profile, add one extra count to avoid 0 count"""
    counting_hash = khmer.new_counting_hash(4, 2000, 1)
    counting_hash.consume(seq)
    composition = [counting_hash.get(kmer)+1 for kmer in kmers]
    if norm == True:
        total = sum(composition)
        composition_norm = [str(number*1.0/total) for number in composition]
        composition = composition_norm
    return composition
    
    
# Functions
def flatten(listOfLists):
    """Flatten a list to one level of nesting"""
    return chain.from_iterable(listOfLists)

def parsemod(dir):
    """Finds the Genemark model, parses it and returns a feature vector with the taxid as the first element"""
    for file in os.listdir(dir):
        if file.endswith("hmm.mod"):
            file = os.path.join(dir,file)
            with open(file, "r") as f:
                rawvec = []
                for line in f.readlines():
                    rawvec.append(line.split())
                start1 = rawvec.index(['COD1']) +1
                stop1 = rawvec.index(['COD2']) -1
                COD1 = (list(flatten(rawvec[start1: stop1 ])))
                start2 = 1 + rawvec.index(['NONC'])
                NONC = (list(flatten(rawvec[start2: start2 +64 ])))
                itemvect = COD1 + NONC
                if len(itemvect) == 256:
                    return itemvect
            f.close()
            break



def main():

    parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
    using emmission data from Metamark')
    parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
    parser.add_argument('--outfile', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
    parser.add_argument('--taxid', help="The taxonomy id")
    parser.add_argument('--label', help="Choice of label, normally taxid, but readid for bining applications", choices=['taxid','readid'],default='taxid')
    parser.add_argument('--mmp', help="the parameters file for GeneMark", default = "../gm_parameters/par_11.modified")
    parser.add_argument('--meta_mmp', help="the parameters file for MetaGeneMark", default = "~/bin/genemark_suite_linux_64/gmsuite/MetaGeneMark_v1.mod")
    parser.add_argument('--tmp', help="root directory to write temp files in", default = "/scratch")
    parser.add_argument('--minlen', help="minimum length to attempt to classify", default = 3000)
    parser.add_argument('--prog', help="GeneMark program to run (all GeneMarkS, GeneMarkS+MetaGenemark, all MetaGeneMark)", choices=['genemarks','metagenemark','hybrid'],default='metagenemark')
    parser.add_argument('--failseq', help="output sequences that failed the program to this file", default = "seq_fail.fa")
    args = parser.parse_args()

    ## File parsing and variable assignment

    mmp = os.path.abspath(args.mmp)
    meta_mmp = os.path.abspath(args.meta_mmp)
    tmp = os.path.abspath(args.tmp)
    file_fail = open(args.failseq,'w')
    
    if not os.path.isdir(tmp): 
        os.mkdir(tmp)
    # for each sequence in the fasta file:
    records = SeqIO.parse(args.input, "fasta")
    cnt_success = 0
    cnt_vectfailure = 0
    cnt_mmfailure = 0
    cnt_probuild_failure = 0
    cnt_gmhmmp_failure = 0
    len_records =0
    shortreads = 0
    kmers = iterate_kmer(4)
    fail_seq = []
    for record in records:
        #go on if reads are too short
        if len(record) < args.minlen:
            shortreads += 1
            continue
        len_records += 1

        tmpdir = tempfile.mkdtemp(dir=tmp) # 
        os.chdir(tmpdir) 
        handle = open("fragment.fasta", "w") # open a fasta file
        SeqIO.write(record, handle, "fasta") # write the sequence to it
        handle.close() # close the file
        if args.prog == "genemarks":
            ## Run GeneMarkS only
            GeneMarkS_params = ["gmsn.pl", "--clean", "--par", mmp,"fragment.fasta"]
            p1 = subprocess.Popen(GeneMarkS_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()
            if p1.returncode == 0: # if genemarkS can generate the model
                featurevect = parsemod(tmpdir)
                if featurevect:
                    if args.label == 'taxid':
                        vect = [args.taxid] + [record.description] + featurevect + \
                        get_composition(str(record.seq).upper(), kmers,True)
                    elif args.label == 'readid':
                        vect = [readid] + [record.description] + featurevect + \
                        get_composition(str(record.seq).upper(), kmers,True)
                    else:
                        raise InputError("the label parameter must be either 'taxid' or 'readid'")
                    args.outfile.write("\t".join(vect))
                    args.outfile.write("\n")
                    shutil.rmtree(tmpdir)
                    cnt_success += 1
                else:
                    #shutil.rmtree(tmpdir)
                    cnt_vectfailure  += 1
                    fail_seq.append(record)
            else: # if not, try to use metagenemark firstly to identify coding region
            
                #shutil.rmtree(tmpdir)
                cnt_mmfailure += 1
                fail_seq.append(record)
                
        elif args.prog == "metagenemark":
            ## Run MetaGenemark only
            MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta"]
            p1 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()
            if p1.returncode == 0: # if MetaGeneMark can generate the gene prediction
            
                probuild_params = ["probuild", "--par", mmp, "--ORDM", "2", "--order_non",\
                 "2", "--revcomp_non", "1", "--seq", "fragment.fasta", "--geneset", \
                 "fragment.fasta.lst", "--mkmod", "hmm.mod"]
                p2 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if probuild can generate hmm model
                
                    featurevect = parsemod(tmpdir)
                    if featurevect:
                        if args.label == 'taxid':
                            vect = [args.taxid] + [record.description] + featurevect + \
                            get_composition(str(record.seq).upper(), kmers,True)
                        elif args.label == 'readid':
                            vect = [readid] + [record.description] + featurevect + \
                            get_composition(str(record.seq).upper(), kmers,True)
                        else:
                            raise InputError("the label parameter must be either 'taxid' or 'readid'")
                        args.outfile.write("\t".join(vect))
                        args.outfile.write("\n")
                        shutil.rmtree(tmpdir)
                        cnt_success += 1
                
                    else:
                        #shutil.rmtree(tmpdir)
                        cnt_vectfailure  += 1
                        fail_seq.append(record)
                else:
                    cnt_probuild_failure += 1
                    fail_seq.append(record)
                    
            else: # if not, try to use metagenemark firstly to identify coding region
            
                #shutil.rmtree(tmpdir)
                cnt_gmhmmp_failure += 1
                fail_seq.append(record)
        elif args.prog == "hybrid":
            ## Run geneMarkS firstly, if failed, (gibbs failed), try metaGenemark, it is still possible
            # that there is no complete cds region to generate hmm model
            GeneMarkS_params = ["gmsn.pl", "--clean", "--par", mmp,"fragment.fasta"]
            p1 = subprocess.Popen(GeneMarkS_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()
            if p1.returncode == 0: # if genemarkS can generate the model
                featurevect = parsemod(tmpdir)
                if featurevect:
                    if args.label == 'taxid':
                        vect = [args.taxid] + [record.description] + featurevect + \
                        get_composition(str(record.seq).upper(), kmers,True)
                    elif args.label == 'readid':
                        vect = [readid] + [record.description] + featurevect + \
                        get_composition(str(record.seq).upper(), kmers,True)
                    else:
                        raise InputError("the label parameter must be either 'taxid' or 'readid'")
                    args.outfile.write("\t".join(vect))
                    args.outfile.write("\n")
                    shutil.rmtree(tmpdir)
                    cnt_success += 1
                else:
                    #shutil.rmtree(tmpdir)
                    cnt_vectfailure  += 1
                    fail_seq.append(record)
            else: # if not, try to use metagenemark firstly to identify coding region
            
                #shutil.rmtree(tmpdir)
                cnt_mmfailure += 1
                
                MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta"]
                p2 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if MetaGeneMark can generate the gene prediction
            
                    probuild_params = ["probuild", "--par", mmp, "--ORDM", "2", "--order_non",\
                     "2", "--revcomp_non", "1", "--seq", "fragment.fasta", "--geneset", \
                     "fragment.fasta.lst", "--mkmod", "hmm.mod"]
                    p3 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                    metamarkout, metamarkerr= p3.communicate()
                    if p3.returncode == 0: # if probuild can generate hmm model
                
                        featurevect = parsemod(tmpdir)
                        if featurevect:
                            if args.label == 'taxid':
                                vect = [args.taxid] + [record.description] + featurevect + \
                                get_composition(str(record.seq).upper(), kmers,True)
                            elif args.label == 'readid':
                                vect = [readid] + [record.description] + featurevect + \
                                get_composition(str(record.seq).upper(), kmers,True)
                            else:
                                raise InputError("the label parameter must be either 'taxid' or 'readid'")
                            args.outfile.write("\t".join(vect))
                            args.outfile.write("\n")
                            shutil.rmtree(tmpdir)
                            cnt_success += 1
                
                        else:
                            #shutil.rmtree(tmpdir)
                            cnt_vectfailure  += 1
                            fail_seq.append(record)
                    else:
                        cnt_probuild_failure += 1
                        fail_seq.append(record)
                    
                else: # 
            
                    #shutil.rmtree(tmpdir)
                    cnt_gmhmmp_failure += 1
                    fail_seq.append(record)
                
    
    SeqIO.write(fail_seq, file_fail, "fasta")
    
    if cnt_success == 0:
        args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, GeneMarkS errors: %s, vector errors: %s, MetaGenemark errors: %s, Probuild errors: %s, reads below min length: %s \n" \
        % (args.taxid, len_records, cnt_success, cnt_mmfailure, cnt_vectfailure, cnt_gmhmmp_failure, cnt_probuild_failure, shortreads))
    args.outfile.close()
    
if __name__ == '__main__':
    main()
