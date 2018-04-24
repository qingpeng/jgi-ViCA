#!/usr/bin/env nextflow
// this script is specifically used for run generating vectors from segments files (large one)
// run locally or on HPC
// if run locally, it will use all the CPUs available to speed up the calculation
// if run on HPC, you can modify nextflow.config to adjust the resource request for HPC

params.segment_file = "/global/projectb/scratch/qpzhang/TARA/Libsvm/all_2k_contigs.fa.tail1000"
params.out = "/global/projectb/scratch/qpzhang/TARA/Libsvm/all_2k_contigs.fa.pfam.tail1000"

params.chunkSize = 6




Channel.fromPath(params.segment_file)
    .splitFasta(by: params.chunkSize, file: true)
    .set { fasta }

process run_pfam {

    input:
    file inputfile from  fasta

    output:
    file 'pfam.out' into vectors

    """
    /global/homes/q/qpzhang/bin/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan --tblout pfam.out /global/homes/q/qpzhang/Pfam_DB/Pfam-A.hmm  ${inputfile}
    """
}

vectors.collectFile(name: params.out)

