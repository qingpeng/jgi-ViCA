#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)
// run locally or on HPC
// if run locally, it will use all the CPUs available to speed up the calculation
// if run on HPC, you can modify nextflow.config to adjust the resource request for HPC

params.out = "./test.fa.14.vect"
params.segment_file = "./test.fa.14"

params.chunkSize = 2

params.package_path = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts"
package_path = Channel.value(params.package_path)

params.genemark_path = "/Users/qingpeng/bin/genemark_suite_macosx/gmsuite/"
genemark_path = Channel.value(params.genemark_path)

params.hmmer_path = "/Users/qingpeng/bin/hmmer-3.1b2-macosx-intel/"
hmmer_path = Channel.value(params.hmmer_path)

params.hmmer_db = "/Users/qingpeng/Pfam_DB/"
hmmer_db = Channel.value(params.hmmer_db)

Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }

process get_feature {
    
    input:
    file inputfile from  fasta
    val package_path
    val genemark_path
    val hmmer_path
    val hmmer_db

    output:
    file 'vector.out' into vectors

    """
    python ${package_path}/6_feature_extraction.py ${inputfile} vector.out ${genemark_path} ${hmmer_path} ${hmmer_db}
    """
}

vectors.collectFile(name: params.out)









