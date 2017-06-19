#!/usr/bin/env python

import argparse
import os
import tempfile
import subprocess
import shutil


def run(genelearn_path, inputfile, outputfile, genemark_path, feature_file,
        model_path, scaler_path):
    tmpdir = os.path.abspath(tempfile.mkdtemp(dir="./"))
    os.chdir(tmpdir)

    vector_file = os.path.basename(inputfile)+'.vect'

    feature_extraction_command = [
        'python',
        genelearn_path+'/6_feature_extraction_noPfam.py', inputfile,
        vector_file, genemark_path]
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

    sklearn_prediction_command = [
        'python',
        genelearn_path + '/sklearn_prediction.py',
        libsvm_file, model_path, scaler_path, 'True', outputfile]

    print "sklearn prediction running...\n"
    return_code = subprocess.call(sklearn_prediction_command)

    if return_code != 0:
        return return_code, "sklearn_prediction"
    print "done!\n"
    shutil.rmtree(tmpdir)
    return 0


def main():

    parser = argparse.ArgumentParser(
        description='A script to do prediction on fasta file')
    parser.add_argument("input_file", help="A multi-sequence fasta file")
    parser.add_argument('output_file', help="output prediction score file")
    parser.add_argument('genemark_path', help="genemark path, ./gmsuite")
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
    feature_file = args.feature_file
    model_path = args.model_directory
    scaler_path = args.scaler_directory

    return_code = run(genelearn_path, input_file, output_file, genemark_path,
                      feature_file, model_path, scaler_path)

    if return_code == 0:
        print "Done!"
    else:
        print return_code


if __name__ == '__main__':
    main()
