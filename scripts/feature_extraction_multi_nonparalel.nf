#!/usr/bin/env nextflow


params.query = "test.fasta"
params.out = "./test.fasta.vect"
params.chunkSize = 10
 
 
Channel
    .fromPath(params.query)
    .splitFasta(by: params.chunkSize)
    .into { fasta }
 


process get_feature {
    
    input:
    file inputfile from  fasta

    output:
    file 'vector.out' into vectors


    """
    feature_genemark_pfam_vfam.py --input ${inputfile} --output_prefix ${inputfile}
    feature_kmer.py  --input ${inputfile} --output ${inputfile}.kmer --ksize 4
    feature_combine.py --output vector.out --length 1  ${inputfile}.kmer ${inputfile}.genemark ${inputfile}.pfam ${inputfile}.vfam
    """
}

vectors.collectFile(name: params.out)



