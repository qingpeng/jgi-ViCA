#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)
// run locally or on HPC


params.out = "./test.fa.vect"
params.segment_file = "./test.fa"

params.chunkSize = 10 
// params.config = "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/config.json.template"
params.package_path = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts"
package_path = Channel.value(params.package_path)


Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }



process get_feature {

    //errorStrategy 'ignore'
    
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

vectors.collectFile(name: params.out)









