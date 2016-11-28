#!/usr/bin/env python


# feature:
# 3.  get lowest e-value for pfam and vfam, if there are multiple hits, move such functions from convert to libsvm to here (09/23 updated)
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
                
                    return zip(range(256),itemvect)
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


def parse_hmmer(file): # extract vector numbers from hmmer result file
# also pick the lowest e value if there are multiple hits
    rawvec = []
    number_processed = 0

    file_obj = open(file,'r')

    e_values = {}
    for line in file_obj:
        if line[0] != "#":
            line = line.rstrip()
            fields = line.split()
            label = fields[0]
            
#            labels.append(fields[0]) # target name
#            values.append(fields[4])
            if not label in e_values:
                e_values[label] = float(fields[4])
            else:
                if float(fields[4])<e_values[label]:
                    e_values[label] = float(fields[4])
    
    
    file_obj.close()
    return e_values.items()


def generate_line(zip_list):
    """ generate line to printout from list of tuples with feature_name and value"""
    line = ''
    for tuple in zip_list:
        pair = str(tuple[0])+":"+ str(tuple[1])
        line = line + ' ' + pair
    return line[1:]
    


def main():

    parser = argparse.ArgumentParser(description='A script to generate feature vector of genemark, pfam and vfam optionally')
    parser.add_argument('--input', help="A multi-sequence fasta file")
    parser.add_argument('--output_prefix', help= "prefix of output file, space delimited format, *.genemark, *.pfam, *.vfam")
    parser.add_argument('--tmp', help="tmp directory", default = "./")
    parser.add_argument('--meta_mmp', help="the parameters file for MetaGeneMark", default = "/global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite/MetaGeneMark_v1.mod")
    parser.add_argument('--mmp', help="the parameters file for GeneMark", default = "/global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/gm_parameters/par_11.modified")
    parser.add_argument('--pfam_db', help="the pfam database, *.hmm", default = "/global/homes/q/qpzhang/Pfam_DB/Pfam-A.hmm")
    parser.add_argument('--vfam_db', help="the vfam database, *.hmm", default = "/global/homes/q/qpzhang/Pfam_DB/vFam-B_2014.hmm")
    parser.add_argument('--IMG_db', help="the IMG virus database, *.hmm", default = "/global/homes/q/qpzhang/Pfam_DB/final_list.hmms")
    parser.add_argument('--feature', help="feature to calculate, (all, genemark) default: all", default = "all")

# Output format:
# seq_id'\t'seq_length'\t'seq_des'\t'vectors, separated by " "
#

    args = parser.parse_args()

    ## File parsing and variable assignment
    mmp = os.path.abspath(args.mmp)
    meta_mmp = os.path.abspath(args.meta_mmp)
    tmp = os.path.abspath(args.tmp)
    pfam_db = os.path.abspath(args.pfam_db)
    vfam_db = os.path.abspath(args.vfam_db)
    img_db = os.path.abspath(args.IMG_db)
    records =  SeqIO.parse(args.input, "fasta")
   # print len(records)
 #   print "here"
    
    output_prefix = args.output_prefix
    
    if args.feature == 'all':
        output_genemark_obj = open(output_prefix+'.genemark','w')
        output_pfam_obj = open(output_prefix+'.pfam','w')
        output_vfam_obj = open(output_prefix+'.vfam','w')
        output_img_obj = open(output_prefix+'.img','w')
    elif args.feature == 'genemark':
        output_genemark_obj = open(output_prefix+'.genemark','w')
        
        

    
    for record in records:
    

        tmpdir = tempfile.mkdtemp(dir=tmp) # 
        os.chdir(tmpdir) 
        handle = open("fragment.fasta", "w") # open a fasta file
        SeqIO.write(record, handle, "fasta") # write the sequence to it
        handle.close() # close the file
        
        length = len(record.seq)
        
  #      print '1'
        if record.description == '':
            des = 'DESCRIPTION'
        else:
            des = record.description
            
        line_genemark = record.id+'\t'+str(length)+'\t'+des+'\t'
        line_pfam = record.id+'\t'+str(length)+'\t'+des+'\t'
        line_vfam = record.id+'\t'+str(length)+'\t'+des+'\t'
        line_img = record.id+'\t'+str(length)+'\t'+des+'\t'
        MetaGeneMark_params = ["gmhmmp", "-m", meta_mmp, "fragment.fasta", "-a", "-A", "fragment.fasta.aa"]
        p1 = subprocess.Popen(MetaGeneMark_params, stdout=subprocess.PIPE)
        metamarkout, metamarkerr= p1.communicate()
        
        if p1.returncode == 0: # if MetaGeneMark can generate the gene prediction 
                probuild_params = ["probuild", "--par", mmp, "--ORDM", "2", "--order_non",\
                 "2", "--revcomp_non", "1", "--seq", "fragment.fasta", "--geneset", \
                 "fragment.fasta.lst", "--mkmod", "hmm.mod"]
                p2 = subprocess.Popen(probuild_params, stdout=subprocess.PIPE)
                metamarkout, metamarkerr= p2.communicate()
                if p2.returncode == 0: # if probuild can generate hmm model
 #                   print "probuild ok\n"
                    vect_genemark = parsemod(tmpdir)
                    line_genemark = line_genemark + generate_line(vect_genemark)
                    
                    if args.feature == 'all':
                    
                        if vect_genemark:
                    
                            hmmscan_params_pfam = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan_pfam","-E",\
                             "1e-5", pfam_db, "fragment.fasta.aa"]
                            p3 = subprocess.Popen(hmmscan_params_pfam, stdout=subprocess.PIPE)
                            metamarkout, metamarkerr= p3.communicate()
                        
                            if p3.returncode == 0:
                                vector_pfam = parse_hmmer("fragment.fasta.aa.hmmscan_pfam")
                                line_pfam = line_pfam + generate_line(vector_pfam)


                            hmmscan_params_vfam = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan_vfam","-E",\
                             "1e-5", vfam_db, "fragment.fasta.aa"]
                            p4 = subprocess.Popen(hmmscan_params_vfam, stdout=subprocess.PIPE)
                            metamarkout, metamarkerr= p4.communicate()
                        
                        
                            if p4.returncode == 0:
                                vector_vfam = parse_hmmer("fragment.fasta.aa.hmmscan_vfam")
                                line_vfam = line_vfam + generate_line(vector_vfam)

                            hmmscan_params_img = ["hmmscan", "--tblout", "fragment.fasta.aa.hmmscan_img","-E",\
                             "1e-5", img_db, "fragment.fasta.aa"]
                            p5 = subprocess.Popen(hmmscan_params_img, stdout=subprocess.PIPE)
                            metamarkout, metamarkerr= p5.communicate()
                        
                        
                            if p5.returncode == 0:
                                vector_img = parse_hmmer("fragment.fasta.aa.hmmscan_img")
                                line_img = line_img + generate_line(vector_img)
            
        output_genemark_obj.write(line_genemark+'\n')
        if args.feature == 'all':
            output_pfam_obj.write(line_pfam+'\n')
            output_vfam_obj.write(line_vfam+'\n')
            output_img_obj.write(line_img+'\n')
        shutil.rmtree(tmpdir)
    
    output_genemark_obj.close()
    if args.feature == 'all':
        output_pfam_obj.close()
        output_vfam_obj.close()
        output_img_obj.close()
    
if __name__ == '__main__':
    main()