#!/bin/sh
#echo "Running shell script"
GENELEARN_PATH="/global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts"

while read LINE; do
    echo ${LINE}
    echo ">test" > segment.fa
    echo ${LINE} >> segment.fa

    python ${GENELEARN_PATH}/2_feature_genemark_pfam_vfam.py --input segment.fa --output_prefix segment.fa
    python ${GENELEARN_PATH}/3_feature_kmer.py  --input segment.fa --output segment.fa.kmer --ksize 4
    python ${GENELEARN_PATH}/4_feature_combine.py --output vector.out --length 1  segment.fa.kmer segment.fa.genemark segment.fa.pfam segment.fa.vfam segment.fa.img 

   cat vector.out
done