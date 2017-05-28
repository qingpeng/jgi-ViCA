#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)
// run locally or on HPC
// if run locally, it will use all the CPUs available to speed up the calculation
// if run on HPC, you can modify nextflow.config to adjust the resource request for HPC

params.out = "./head10000.fasta.vect"
params.segment_file = "./head10000.fasta"

params.chunkSize = 10

params.package_path = "/global/homes/q/qpzhang/Github/jgi-ViCA/scripts"
package_path = Channel.value(params.package_path)

params.genemark_path = "/global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite"
genemark_path = Channel.value(params.genemark_path)



Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }

process get_feature {
    
    input:
    file inputfile from  fasta
    val package_path
    val genemark_path

    output:
    file 'vector.out' into vectors

    """
    python ${package_path}/6_feature_extraction_noPfam.py ${inputfile} vector.out ${genemark_path}
    """
}

vectors.collectFile(name: params.out)









