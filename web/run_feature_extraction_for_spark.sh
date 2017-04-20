#!/bin/sh
#echo "Running shell script"
GENELEARN_PATH = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts/"
GENEMARK_PATH = "~/bin/genemark_suite_macosx/gmsuite/"
HMMER_PATH = "~/bin/hmmer-3.1b2-macosx-intel/"
PFAM_DB = "~/Local/Pfam_DB/"


source /Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts/Flask/env.sh

while read LINE; do
    #echo ${LINE}
    
    echo ">test" > segment.fa
    #echo $HOME > env.fa
    echo ${LINE} >> segment.fa

    python ${GENELEARN_PATH}/2_feature_genemark_pfam_vfam.py --input segment.fa --output_prefix segment.fa
    python ${GENELEARN_PATH}/3_feature_kmer.py  --input segment.fa --output segment.fa.kmer --ksize 4
    python ${GENELEARN_PATH}/4_feature_combine.py --output vector.out --length 1  segment.fa.kmer segment.fa.genemark segment.fa.pfam segment.fa.vfam segment.fa.img 
    python ${GENELEARN_PATH}/vect2svmlib.py -v vector.out -f ${GENELEARN_PATH}/nonvirus.vect.feature_index.new_list -o vector.out.smvlib -t 1
    python ${GENELEARN_PATH}/modify_label.py vector.out.smvlib 2104 >vector.out.smvlib2

python ${GENELEARN_PATH}/6_feature_extraction.py virus_nonvirus_3seqs.fa virus_nonvirus_3seqs.fa.vect ~/bin/genemark_suite_macosx/gmsuite/ ~/bin/hmmer-3.1b2-macosx-intel/ ~/Local/Pfam_DB/


python ${GENELEARN_PATH}/prepare_libsvm_for_prediction.py virus_nonvirus_3seqs.fa.vect ~/Dropbox/Development/Github/jgi-ViCA/scripts/model/all_segment.fasta.vect.feature_index virus_nonvirus_3seqs.fa.vect.libsvm

    cat vector.out.smvlib2
done