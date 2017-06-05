#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)
// run locally or on HPC
// if run locally, it will use all the CPUs available to speed up the calculation
// if run on HPC, you can modify nextflow.config to adjust the resource request for HPC

libsvm_files = Channel.fromPath('/global/projectb/scratch/qpzhang/TARA/Libsvm/*.libsvm')




Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }

process get_feature {
    
    input:
    file libsvm from  libsvm_files

    """
    python /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/sklearn_prediction.py ${libsvm} /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/model/Sklearn/all_segment.fasta.vect.family.training.svmlib.no4.0_1_scikit_model_median_True_1x /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/model/Sklearn/all_segment.fasta.vect.family.training.svmlib.no4.0_1_scikit_scaler_median_True_1x True ${libsvm}.prediction_sklearn
    """
}








