#!/usr/bin/env python

# v3, add support to use pfam to generate vectors 03/04/2016
# 

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

# test

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
        
        
def get_composition(ksize, seq, kmers, norm):
    """ get the composition profile, add one extra count to avoid 0 count"""
    counting_hash = khmer.new_counting_hash(ksize, 200000, 1)
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
# 
# #                                                                               --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# # target name        accession  query name                           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc descript
# #------------------- ----------                 -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- --------
# COX2_TM              PF02790.12 gene_4|GeneMark.hmm|92_aa|+|4616|4894 -            1.5e-15   56.9   0.7   1.9e-15   56.6   0.7   1.1   1   0   0   1   1   1   1 Cytochr
# #
# # Program:         hmmscan
# # Version:         3.1b2 (February 2015)
# #


def parsemod_pfam(dir): # output the hit accession label and the score
    rawvec = []
    number_processed = 0
    for file in os.listdir(dir):
        if file.endswith("fragment.fasta.aa.hmmscan_pfam") or file.endswith("fragment.fasta.aa.hmmscan_vfam"):
            number_processed += 1
            file = os.path.join(dir,file)
            file_obj = open(file,'r')
            
            for line in file_obj:
                if line[0] != "#":
                    line = line.rstrip()
                    fields = line.split()
                    label = fields[0] # target name
                    e_value = fields[4]
                    rawvec.append(label+":"+e_value)
            file_obj.close()
            if number_processed == 2: # found and processed pfam and vfam searching results
                return rawvec
                break
            
    return rawvec
    
            




def main():

    parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
    using emmission data from Metamark')
    parser.add_argument('--input', help="A multi-sequence fasta file",type=argparse.FileType('r'), default='-')
    parser.add_argument('--outfile', help= "Output file, tab delimited format", type=argparse.FileType('w'), default='-')
    parser.add_argument('--taxid', help="The taxonomy id")
    parser.add_argument('--label', help="Choice of label, normally taxid, but readid for bining applications", choices=['taxid','readid'],default='taxid')
    parser.add_argument('--mmp', help="the parameters file for GeneMark", default = "../gm_parameters/par_11.modified")
    parser.add_argument('--meta_mmp', help="the parameters file for MetaGeneMark", default = "/global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite/MetaGeneMark_v1.mod")
    parser.add_argument('--tmp', help="root directory to write temp files in", default = "/scratch")
    parser.add_argument('--minlen', help="minimum length to attempt to classify", default = 1000)
    parser.add_argument('--prog', help="GeneMark program to run ( genemarks - all GeneMarkS, hybrid - GeneMarkS+MetaGenemark, metagenemark - all MetaGeneMark)", choices=['genemarks','metagenemark','hybrid','pfam','pfam_combine','pfam_combine_vfam'],default='metagenemark')
    parser.add_argument('--failseq', help="output sequences that failed the program to this file", default = "seq_fail.fa")
    parser.add_argument('--ksize', help="size of kmer, default=4", default = 4)
    parser.add_argument('--fam_path', help="path with pfam/vfam database files", default = 4)
 
    args = parser.parse_args()

    ## File parsing and variable assignment

    mmp = os.path.abspath(args.mmp)
    meta_mmp = os.path.abspath(args.meta_mmp)
    fam_path = args.fam_path
    tmp = os.path.abspath(args.tmp)
#    file_fail = open(args.failseq,'w')
    
    if not os.path.isdir(tmp): 
        os.mkdir(tmp)
    # for each sequence in the fasta file:
    records = SeqIO.parse(args.input, "fasta")
    cnt_success = 0
    cnt_vectfailure = 0
    cnt_mmfailure = 0
    cnt_probuild_failure = 0
    cnt_gmhmmp_failure = 0
    cnt_hmmscan_failure = 0
    len_records =0
    shortreads = 0
    k_size = int(args.ksize)
    kmers = iterate_kmer(k_size)
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
                        get_composition(k_size,str(record.seq).upper(), kmers,True)
                    elif args.label == 'readid':
                        vect = [readid] + [record.description] + featurevect + \
                        get_composition(k_size,str(record.seq).upper(), kmers,True)
                    else:
                        raise InputError("the label parameter must be either 'taxid' or 'readid'")
                    args.outfile.write("\t".join(vect))
                    args.outfile.write("\n")
                    shutil.rmtree(tmpdir)
                    cnt_success += 1
                else:
                    shutil.rmtree(tmpdir)
                    cnt_vectfailure  += 1
                    fail_seq.append(record)
            else: 
            
                shutil.rmtree(tmpdir)
                cnt_mmfailure += 1
                fail_seq.append(record)
                
        elif args.prog == "metagenemark":
            ## Run MetaGenemark only
            MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta"]
            p1 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()
            if p1.returncode == 0: # if MetaGeneMark can generate the gene prediction
#                print "gmhmmp ok\n"
                probuild_params = ["probuild", "--par", mmp, "--ORDM", "2", "--order_non",\
                 "2", "--revcomp_non", "1", "--seq", "fragment.fasta", "--geneset", \
                 "fragment.fasta.lst", "--mkmod", "hmm.mod"]
                p2 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if probuild can generate hmm model
 #                   print "probuild ok\n"
                    featurevect = parsemod(tmpdir)
                    if featurevect:
                        if args.label == 'taxid':
                            vect = [args.taxid] + [record.description] + featurevect + \
                            get_composition(k_size,str(record.seq).upper(), kmers,True)
                        elif args.label == 'readid':
                            vect = [readid] + [record.description] + featurevect + \
                            get_composition(k_size,str(record.seq).upper(), kmers,True)
                        else:
                            raise InputError("the label parameter must be either 'taxid' or 'readid'")
                        args.outfile.write("\t".join(vect))
                        args.outfile.write("\n")
                        shutil.rmtree(tmpdir)
                        cnt_success += 1
                
                    else:
                        shutil.rmtree(tmpdir)
                        cnt_vectfailure  += 1
                        fail_seq.append(record)
                else:
  #                  print "probuild fail\n"
                    shutil.rmtree(tmpdir)
                    cnt_probuild_failure += 1
                    fail_seq.append(record)
                    
            else: # 
   #             print "gmhmmp fail\n"
                shutil.rmtree(tmpdir)
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
                        get_composition(k_size,str(record.seq).upper(), kmers,True)
                    elif args.label == 'readid':
                        vect = [readid] + [record.description] + featurevect + \
                        get_composition(k_size,str(record.seq).upper(), kmers,True)
                    else:
                        raise InputError("the label parameter must be either 'taxid' or 'readid'")
                    args.outfile.write("\t".join(vect))
                    args.outfile.write("\n")
                    shutil.rmtree(tmpdir)
                    cnt_success += 1
                else:
                    shutil.rmtree(tmpdir)
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
                                get_composition(k_size,str(record.seq).upper(), kmers,True)
                            elif args.label == 'readid':
                                vect = [readid] + [record.description] + featurevect + \
                                get_composition(k_size,str(record.seq).upper(), kmers,True)
                            else:
                                raise InputError("the label parameter must be either 'taxid' or 'readid'")
                            args.outfile.write("\t".join(vect))
                            args.outfile.write("\n")
                            shutil.rmtree(tmpdir)
                            cnt_success += 1
                
                        else:
                            shutil.rmtree(tmpdir)
                            cnt_vectfailure  += 1
                            fail_seq.append(record)
                    else:
                        cnt_probuild_failure += 1
                        fail_seq.append(record)
                        shutil.rmtree(tmpdir)
                    
                else: # 
            
                    shutil.rmtree(tmpdir)
                    cnt_gmhmmp_failure += 1
                    fail_seq.append(record)
          
        elif args.prog == "pfam":

            ## Run MetaGenemark only
            MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta", "-a", "-A", "fragment.fasta.aa"]
            p1 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()
            acc_file_path = "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam-A.hmm.acc"
            pfam_path = "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam-A.hmm"
            acc_file = open(acc_file_path,'r')
            acc_list = []
            for line in acc_file:
                line = line.rstrip()
                acc_list.append(line.split()[1])
            
                
            if p1.returncode == 0: # if MetaGeneMark can generate the gene prediction
#                print "gmhmmp ok\n"

                hmmscan_params = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan","-E",\
                 "1e-5", pfam_path, "fragment.fasta.aa"]
                 

                p2 = subprocess.Popen(hmmscan_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if probuild can generate hmm model
 #                   print "probuild ok\n"
                    
                    featurevect = parsemod_pfam(tmpdir, acc_list)
                    
                    
                    if featurevect:
                        if args.label == 'taxid':
                            vect = [args.taxid] + [record.description] + featurevect
                        elif args.label == 'readid':
                            vect = [readid] + [record.description] + featurevect 
                        else:
                            raise InputError("the label parameter must be either 'taxid' or 'readid'")
                        args.outfile.write("\t".join(vect))
                        args.outfile.write("\n")
                        shutil.rmtree(tmpdir)
                        cnt_success += 1
                
                    else:
                        shutil.rmtree(tmpdir)
                        cnt_vectfailure  += 1
                        fail_seq.append(record)
                else:
  #                  print "probuild fail\n"
                    shutil.rmtree(tmpdir)
                    cnt_hmmscan_failure += 1
                    fail_seq.append(record)
                    
            else: # 
   #             print "gmhmmp fail\n"
                shutil.rmtree(tmpdir)
                cnt_gmhmmp_failure += 1
                fail_seq.append(record)
                
                
                
                
        elif args.prog == "pfam_combine":

            ## Run MetaGenemark only
            MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta", "-a", "-A", "fragment.fasta.aa"]
            p1 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()
            acc_file_path = "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam-A.hmm.acc"
            pfam_path = "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam-A.hmm"
            acc_file = open(acc_file_path,'r')
            acc_list = []
            for line in acc_file:
                line = line.rstrip()
                acc_list.append(line.split()[1])
            
                
            if p1.returncode == 0: # if MetaGeneMark can generate the gene prediction
#                print "gmhmmp ok\n"
                probuild_params = ["probuild", "--par", mmp, "--ORDM", "2", "--order_non",\
                 "2", "--revcomp_non", "1", "--seq", "fragment.fasta", "--geneset", \
                 "fragment.fasta.lst", "--mkmod", "hmm.mod"]
                p2 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if probuild can generate hmm model
 #                   print "probuild ok\n"
                    featurevect_probuild = parsemod(tmpdir)
                    if featurevect_probuild:
                    

                        hmmscan_params = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan","-E",\
                         "1e-5", pfam_path, "fragment.fasta.aa"]
                 

                        p3 = subprocess.Popen(hmmscan_params, stdout=subprocess.PIPE)
                        metamarkout, metamarkerr= p3.communicate()
                        if p3.returncode == 0: # if hmmscan can generate results
         #                   print "probuild ok\n"
                    
                            featurevect_hmm = parsemod_pfam(tmpdir, acc_list)
                            
                            if featurevect_hmm:
                                featurevect_hmm_out = featurevect_hmm
                            else:
                                featurevect_hmm_out = []
                                cnt_hmmscan_failure += 1
                                

                        else:
                            featurevect_hmm_out = []
                            cnt_hmmscan_failure += 1
                            
                        if args.label == 'taxid':
                            vect = [args.taxid] + [record.description] + featurevect_probuild + \
                            get_composition(k_size,str(record.seq).upper(), kmers,True) + featurevect_hmm_out
                        elif args.label == 'readid':
                            vect = [readid] + [record.description] + featurevect_probuild + \
                            get_composition(k_size,str(record.seq).upper(), kmers,True) + featurevect_hmm_out
                        else:
                            raise InputError("the label parameter must be either 'taxid' or 'readid'")
                        args.outfile.write("\t".join(vect))
                        args.outfile.write("\n")
                        shutil.rmtree(tmpdir)
                        cnt_success += 1
                                

                
                    else:
                        shutil.rmtree(tmpdir)
                        cnt_vectfailure  += 1
                        fail_seq.append(record)
                else:
  #                  print "probuild fail\n"
                    shutil.rmtree(tmpdir)
                    cnt_probuild_failure += 1
                    fail_seq.append(record)
                    
            else: # 
   #             print "gmhmmp fail\n"
                shutil.rmtree(tmpdir)
                cnt_gmhmmp_failure += 1
                fail_seq.append(record)
                
        elif args.prog == "pfam_combine_vfam":

            ## Run MetaGenemark only
            MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta", "-a", "-A", "fragment.fasta.aa"]
            p1 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
            metamarkout, metamarkerr= p1.communicate()

            #pfam_path = "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Pfam-A.hmm"
            pfam_hmm = fam_path + '/Pfam-A.hmm'
            vfam_hmm = fam_path + '/vFam-B_2014.hmm'
            #vfam_path = "/global/projectb/scratch/qpzhang/Full_Training/Pfam/Vfam/vFam-B_2014.hmm"

            
                
            if p1.returncode == 0: # if MetaGeneMark can generate the gene prediction
#                print "gmhmmp ok\n"
                probuild_params = ["probuild", "--par", mmp, "--ORDM", "2", "--order_non",\
                 "2", "--revcomp_non", "1", "--seq", "fragment.fasta", "--geneset", \
                 "fragment.fasta.lst", "--mkmod", "hmm.mod"]
                p2 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if probuild can generate hmm model
 #                   print "probuild ok\n"
                    featurevect_probuild = parsemod(tmpdir)
                    if featurevect_probuild:
                    

                        hmmscan_params_pfam = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan_pfam","-E",\
                         "1e-5", pfam_hmm, "fragment.fasta.aa"]
                        p3 = subprocess.Popen(hmmscan_params_pfam, stdout=subprocess.PIPE)
                        metamarkout, metamarkerr= p3.communicate()
                        
                        hmmscan_params_vfam = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan_vfam","-E",\
                         "1e-5", vfam_hmm, "fragment.fasta.aa"]
                        p4 = subprocess.Popen(hmmscan_params_vfam, stdout=subprocess.PIPE)
                        metamarkout, metamarkerr= p4.communicate()
                        
                        
                        
                        if p3.returncode == 0 or p4.returncode == 0: # if hmmscan can generate results
                    
                            featurevect_hmm = parsemod_pfam(tmpdir)
                            
                            if featurevect_hmm:
                                featurevect_hmm_out = featurevect_hmm
                            else:
                                featurevect_hmm_out = []
                                cnt_hmmscan_failure += 1
                                

                        else:
                            featurevect_hmm_out = []
                            cnt_hmmscan_failure += 1
                            
                        if args.label == 'taxid':
                            vect = [args.taxid] + [record.description] + featurevect_probuild + \
                            get_composition(k_size,str(record.seq).upper(), kmers,True) + featurevect_hmm_out
                        elif args.label == 'readid':
                            vect = [readid] + [record.description] + featurevect_probuild + \
                            get_composition(k_size,str(record.seq).upper(), kmers,True) + featurevect_hmm_out
                        else:
                            raise InputError("the label parameter must be either 'taxid' or 'readid'")
                        args.outfile.write("\t".join(vect))
                        args.outfile.write("\n")
                        shutil.rmtree(tmpdir)
                        cnt_success += 1
                                

                
                    else:
                        shutil.rmtree(tmpdir)
                        cnt_vectfailure  += 1
                        fail_seq.append(record)
                else:
  #                  print "probuild fail\n"
                    shutil.rmtree(tmpdir)
                    cnt_probuild_failure += 1
                    fail_seq.append(record)
                    
            else: # 
   #             print "gmhmmp fail\n"
                shutil.rmtree(tmpdir)
                cnt_gmhmmp_failure += 1
                fail_seq.append(record)
                

#    SeqIO.write(fail_seq, file_fail, "fasta")
    
    if cnt_success == 0:
        args.outfile.write("#Taxon id: %s, Number of Contigs: %s, Successes: %s, GeneMarkS errors: %s, vector errors: %s, MetaGenemark errors: %s,hmmscan errors: %s,  Probuild errors: %s, reads below min length: %s \n" \
        % (args.taxid, len_records, cnt_success, cnt_mmfailure, cnt_vectfailure, cnt_gmhmmp_failure, cnt_hmmscan_failure, cnt_probuild_failure, shortreads))
    args.outfile.close()
    
if __name__ == '__main__':
    main()
