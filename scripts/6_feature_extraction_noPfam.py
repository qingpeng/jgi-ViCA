#!/usr/bin/env python

import argparse
import os
import tempfile
import subprocess
import shutil


def run(genelearn_path, inputfile, outputfile, genemark_path):
    tmpdir = os.path.abspath(tempfile.mkdtemp(dir="./"))
    os.chdir(tmpdir)

    feature_genemark_command = [
        genelearn_path+'/2_feature_genemark.py',
        '--input', inputfile,
        '--output_prefix', "output", '--genemark_path', genemark_path]
    return_code = subprocess.call(feature_genemark_command)

    if return_code != 0:
        return return_code, 2

    feature_kmer_command = [
        genelearn_path+'/3_feature_kmer.py', '--input', inputfile,
        '--output', 'output.kmer', '--ksize', '4']

    return_code = subprocess.call(feature_kmer_command)

    if return_code != 0:
        return return_code, 3

    feature_combine_command = [
        genelearn_path+'/4_feature_combine.py', '--output', outputfile,
        '--length', '1', 'output.kmer', 'output.genemark']

    return_code = subprocess.call(feature_combine_command)

    if return_code != 0:
        return return_code, 4

    shutil.rmtree(tmpdir)
    return 0


def main():

    parser = argparse.ArgumentParser(
        description='A script to generate feature vector')
    parser.add_argument("input_file", help="A multi-sequence fasta file")
    parser.add_argument('output_file', help="output vector file")
    parser.add_argument('genemark_path', help="genemark path, ./gmsuite")

    args = parser.parse_args()

    genelearn_path = os.path.dirname(os.path.realpath(__file__))

    input_file = os.path.abspath(args.input_file)
    output_file = os.path.abspath(args.output_file)
    genemark_path = args.genemark_path

    return_code = run(genelearn_path, input_file, output_file, genemark_path)

    if return_code == 0:
        print "Done!"
    else:
        print return_code


if __name__ == '__main__':
    main()