#!/usr/bin/env nextflow


params.query = "./Output/{1,2,3,4,5,6,7,8,9,0}*"
params.out = "./output.fasta.vect2"
params.chunkSize = 50
 
 
// run create_shred.py to get shreded segments for the genomes under taxonomy id

// run reftree.py to get the fasta file from the reftree db

//


 
 
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



