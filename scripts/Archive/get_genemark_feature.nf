#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)



params.out = "/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/Split_with_name/training_25k.vect.fasta.genemark_feature"
params.segment_file = "/global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/Split_with_name/training_25k.vect.fasta"

params.chunkSize = 2000

params.package_path = "/global/homes/q/qpzhang/Github/jgi-ViCA/scripts"
package_path = Channel.value(params.package_path)


Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }


process get_feature {

    errorStrategy 'ignore'
    
    input:
    file inputfile from  fasta
    val package_path

    output:
    file ${inputfile}.genemark into vectors

    """
    python ${package_path}/2_feature_genemark_pfam_vfam.py --input ${inputfile} --output_prefix ${inputfile} --genemark_path /global/homes/q/qpzhang/bin/genemark_suite_linux_64/gmsuite

    """
}

vectors.collectFile(name: params.out)









