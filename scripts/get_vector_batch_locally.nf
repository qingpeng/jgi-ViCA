#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)
// run locally or on HPC


params.out = "./test.fa.14.vect"
params.segment_file = "./test.fa.14"

params.chunkSize = 2
// params.config = "/global/projectb/scratch/qpzhang/Run_Genelearn/Small_Set/Nextflow/config.json.template"
params.package_path = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/scripts"
package_path = Channel.value(params.package_path)

params.genemark_path = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/ThirdParty/genemark_suite_macosx/gmsuite/"
genemark_path = Channel.value(params.genemark_path)

params.hmmer_path = "/Users/qingpeng/Dropbox/Development/Bitbucket/jgi-genelearn/ThirdParty/hmmer-3.1b2-macosx-intel/"
hmmer_path = Channel.value(params.hmmer_path)

params.hmmer_db = "/Users/qingpeng/Pfam_DB/"
hmmer_db = Channel.value(params.hmmer_db)


Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }



process get_feature {

    //errorStrategy 'ignore'
    
    input:
    file inputfile from  fasta
    val package_path
    val genemark_path
    val hmmer_path
    val hmmer_db

    output:
    file 'vector.out' into vectors


    """
    python ${package_path}/2_feature_genemark_pfam_vfam.py --input ${inputfile} --output_prefix ${inputfile} --genemark_path ${genemark_path} --hmmer_path ${hmmer_path} --hmmer_db ${hmmer_db}
    python ${package_path}/3_feature_kmer.py  --input ${inputfile} --output ${inputfile}.kmer --ksize 4
    python ${package_path}/4_feature_combine.py --output vector.out --length 1  ${inputfile}.kmer ${inputfile}.genemark ${inputfile}.pfam ${inputfile}.vfam ${inputfile}.img
    """
}

vectors.collectFile(name: params.out)









