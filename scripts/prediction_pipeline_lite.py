#!/usr/bin/env python

import argparse
import os
import tempfile
import subprocess
import shutil


def run(genelearn_path, inputfile, outputfile, genemark_path, hmmer_path,
        hmmer_db, spark_path, feature_file, model_path, scaler_path):
    tmpdir = os.path.abspath(tempfile.mkdtemp(dir="./"))
    os.chdir(tmpdir)

    vector_file = os.path.basename(inputfile)+'.vect'

    feature_extraction_command = [
        'python',
        genelearn_path+'/6_feature_extraction.py', inputfile,
        vector_file, genemark_path, hmmer_path, hmmer_db]
    print feature_extraction_command
    print "feature extraction running...\n"
    return_code = subprocess.call(feature_extraction_command)

    if return_code != 0:
        return return_code, "feature_extraction"
    # exit()

    libsvm_file = vector_file+'.libsvm'
    prepare_libsvm_command = [
        'python',
        genelearn_path+'/prepare_libsvm_for_prediction.py', vector_file,
        feature_file, libsvm_file]
    print "prepare libsvm running...\n"
    return_code = subprocess.call(prepare_libsvm_command)

    if return_code != 0:
        return return_code, "prepare_libsvm"

    libsvm_file_no4 = libsvm_file+'.no4'  # without HMM vfam features
    prepare_libsvm_command = [
        'python',
        genelearn_path+'/pick_vectors_by_feature.py',
        '-d', genelearn_path+'/model/all_segment.fasta.vect.feature_index',
        '-f', '0_1_2_3', '-i', libsvm_file, '-o', libsvm_file_no4]
    print "pick features from libsvm running...\n"
    return_code = subprocess.call(prepare_libsvm_command)

    if return_code != 0:
        return return_code, "pick_features"

    spark_prediction_command = [
        spark_path+'bin/spark-submit',
        genelearn_path+'/spark_prediction_dataframe.py', libsvm_file_no4,
        model_path,
        scaler_path, outputfile]
    print "spark prediction running...\n"
    return_code = subprocess.call(spark_prediction_command)

    if return_code != 0:
        return return_code, "spark_prediction"
    print "done!\n"
    shutil.rmtree(tmpdir)
    return 0


def main():

    parser = argparse.ArgumentParser(
        description='A script to do prediction on fasta file')
    parser.add_argument("input_file", help="A multi-sequence fasta file")
    parser.add_argument('output_file', help="output prediction score file")
    parser.add_argument('genemark_path', help="genemark path, ./gmsuite")
    parser.add_argument('hmmer_path', help="HMMER path")
    parser.add_argument('hmmer_db', help="HMMER db, with pfam, vfam, imgdb")
    parser.add_argument('spark_path',
                        help="Spark path, spark-2.1.0-bin-hadoop2.7")
    parser.add_argument('feature_file',
                        help="feature list file from training data")
    parser.add_argument('model_directory',
                        help="model used for prediction")
    parser.add_argument('scaler_directory',
                        help="scaler used for normalization")
    args = parser.parse_args()

    genelearn_path = os.path.dirname(os.path.realpath(__file__))

    input_file = os.path.abspath(args.input_file)
    output_file = os.path.abspath(args.output_file)
    genemark_path = args.genemark_path
    hmmer_path = args.hmmer_path
    hmmer_db = args.hmmer_db
    spark_path = args.spark_path
    feature_file = args.feature_file
    model_path = args.model_directory
    scaler_path = args.scaler_directory

    return_code = run(genelearn_path, input_file, output_file, genemark_path,
                      hmmer_path, hmmer_db, spark_path, feature_file,
                      model_path, scaler_path)

    if return_code == 0:
        print "Done!"
    else:
        print return_code


if __name__ == '__main__':
    main()
