#!/usr/bin/env nextflow


myFileChannel1 = Channel.fromPath('test.fasta')
myFileChannel2= Channel.fromPath('test.fasta')

process genemark_pfam_vfam {
    input:
    file x_genemark_pfam_vfam from  myFileChannel1

    output:
    file "test.fasta.genemark" into genemark
    file "test.fasta.pfam" into pfam
    file "test.fasta.vfam" into vfam

    """
    feature_genemark_pfam_vfam.py --input ${x_genemark_pfam_vfam} --output_prefix test.fasta
    """
}

process kmer {
    input:
    file x_kmer from  myFileChannel2 

    output:
    file "test.fasta.kmer" into kmer
    
    """
    feature_kmer.py  --input ${x_kmer} --output test.fasta.kmer --ksize 4
    """
}


process feature_combine {
    publishDir './'
    input:
    file in_genemark from genemark
    file in_pfam from pfam
    file in_vfam from vfam
    file in_kmer from kmer
    
    output:
    file "all.vect" into vect
    
    """
    feature_combine.py --output all.vect --length 1 ${in_genemark} ${in_pfam} ${in_vfam} ${in_kmer}
    """
}


