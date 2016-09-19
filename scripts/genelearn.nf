#!/usr/bin/env nextflow



params.out = "./output.fasta.vect_uge2"
params.chunkSize = 50
params.config = "/global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/config.json.template"
 
// run create_shred.py to get shreded segments for the genomes under taxonomy id

// run reftree.py to get the fasta file from the reftree db

//


configfile = Channel.fromPath( params.config )


process get_segment {
    
    input:
    file configfile
    
    output:
    file "Output/{1,2,3,4,5,6,7,8,9,0}*" into segment
    
    """
    create_shred.py -c $configfile -o Output
    """
}

 
segment
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



