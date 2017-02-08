#!/usr/bin/env nextflow



params.out = "./fasta.vect"
params.segment_file = "/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/all_segment.fa"
// params.segment_file = "/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/50000.fa"

params.chunkSize = 2400
// params.config = "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/config.json.template"
params.package_path = "/global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts"
package_path = Channel.value(params.package_path)

// run create_shred.py to get shreded segments for the genomes under taxonomy id

// run reftree.py to get the fasta file from the reftree db

//


Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize)
    .into { fasta }



process get_feature {

    input:
    file inputfile from  fasta
    val package_path

    output:
    file 'vector.out' into vectors


    """
    python ${package_path}/2_feature_genemark_pfam_vfam.py --input ${inputfile} --output_prefix ${inputfile}
    python ${package_path}/3_feature_kmer.py  --input ${inputfile} --output ${inputfile}.kmer --ksize 4
    python ${package_path}/4_feature_combine.py --output vector.out --length 1  ${inputfile}.kmer ${inputfile}.genemark ${inputfile}.pfam ${inputfile}.vfam ${inputfile}.img
    """
}

vectors.collectFile(name: params.out).into {vectors_combine}



process split_vectors_on_rank {
    publishDir "results"

    input:
    file vectors_combine
    val package_path

    output:
    file "${vectors_combine}.*.svmlib" into svmlib
    file "*.feature_index" into feature_index
    // feature_index with the list of feature

    """
    python ${package_path}/5_create_training_testing.py -v ${vectors_combine} -r True
    """
}











