#!/usr/bin/env nextflow
 
/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */
params.query = "test.fa"
params.out = "./test.fa.vect"
params.chunkSize = 10
 
 
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta'
 */
Channel
    .fromPath(params.query)
    .splitFasta(by: params.chunkSize)
    .set { fasta }
 
/*
 * Executes a BLAST job for each chunk emitted by the 'fasta' channel
 * and creates as output a channel named 'top_hits' emitting the resulting
 * BLAST matches 
 */
process blast {
    input:
    file 'query.fa' from fasta
 
    output:
    file 'vect' into vector
 
    """
    feature_extraction_v2.py --input query.fa  --outfile vect --label readid --mmp ~/Bitbucket/jgi-genelearn/scripts/gm_parameters/par_11.modified --prog pfam_combine_vfam --fam_path ~/Pfam_DB/
    """
}
 
 
 
/*
 * Collects all the sequences files into a single file
 * and prints the resulting file content when complete
 */
vector.collectFile(name: params.out)
//    .println { file -> "matching sequences:\n ${file.text}" }

