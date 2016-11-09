#!/usr/bin/env nextflow

// this script is used for testing the pipeline, 
// using a non-virus subtree, 
// split by the rank using new Usearch approach
// using virus node, does not split , just convert from vector into libsvm 
// format, this will be similar to real data to classify
//


params.out = "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/nonvirus.vect"
params.chunkSize = 52
params.config = "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/config.json.template"
params.package_path = "/global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts"

 
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
    python ${params.package_path}/1_create_shred.py -c $configfile -o Output
    """
}


// cat */* >all_segment.fa &


segment
    .flatMap().splitFasta(by: params.chunkSize, file: true)
    .set { fasta }
 


process get_feature {
    
    input:
    file inputfile from  fasta


    output:
    file 'vector.out' into vectors


    """
    python ${params.package_path}/2_feature_genemark_pfam_vfam.py --input ${inputfile} --output_prefix ${inputfile}
    python ${params.package_path}/3_feature_kmer.py  --input ${inputfile} --output ${inputfile}.kmer --ksize 4
    python ${params.package_path}/4_feature_combine.py --output vector.out --length 1  ${inputfile}.kmer ${inputfile}.genemark ${inputfile}.pfam ${inputfile}.vfam ${inputfile}.img 
    """
}

vectors.collectFile(name: params.out).into {vectors_combine}



process split_vectors_on_rank {
    publishDir "nonvirus_results"

    input:
    file vectors_combine

 
    output:
    file "${vectors_combine}.*.svmlib" into svmlib
    file "*.feature_index" into feature_index
    // feature_index with the list of feature 
    
    """
    python ${params.package_path}/5_create_training_testing.py -v ${vectors_combine} -r True
    """
}
    
    
    
    
    
    
    
    
    
    
    
