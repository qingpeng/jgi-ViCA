#!/usr/bin/env nextflow



params.out = "./69657.fasta.vect"
params.chunkSize = 52
params.config = "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/config.json.template"
 
// run create_shred.py to get shreded segments for the genomes under taxonomy id

// run reftree.py to get the fasta file from the reftree db

//


configfile = Channel.fromPath( params.config )


process get_segment {
    
    input:
    file configfile
    
    output:
    file "Output/[0-9]*" into segment
    
    """
    create_shred.py -c $configfile -o Output
    """
}

 
segment
    .flatMap().splitFasta(by: params.chunkSize)
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



