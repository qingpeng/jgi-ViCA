#!/usr/bin/env python

import argparse


class Feature:
    def __init__(self, index_file):
        self.index_file = index_file

    def get_list_to_extract(self, feature_code):
        feature_code_set = set(feature_code.split('_'))
        file_index_file_obj = open(self.index_file, 'r')
        feature_id_set = set()
        for line in file_index_file_obj:
            line = line.rstrip()
            fields = line.split()
            feature_code = fields[1][0]
            if feature_code in feature_code_set:
                feature_id_set.add(int(fields[0])+1)
                # feature list starts from 0
        return feature_id_set


class VectorFile:
    def __init__(self, vector_file):
        self.vector_file = vector_file

    @staticmethod
    def process_line(line, feature_id_set):

        fields = line.split()
        newline = fields[0]
        for vector in fields[1:]:
            items = vector.split(":")
            if int(items[0]) in feature_id_set:
                newline = newline + ' ' + vector
        return newline

    def extract_feature(self, feature_id_set, outputfile):
        file_vector_file_obj = open(self.vector_file, 'r')
        file_output_obj = open(outputfile, 'w')

        for line in file_vector_file_obj:
            line = line.rstrip()

            new_line = self.process_line(line, feature_id_set)

            file_output_obj.write(new_line+'\n')

        file_vector_file_obj.close()
        file_output_obj.close()


def main():
    parser = argparse.ArgumentParser(
        description='A script to pick vectors by feature id')
    parser.add_argument('-d', '--index',
                        help='feature_index file', required=True)
    parser.add_argument(
        '-f', '--feature', help='feature id to pick, format: 0_1_2_3, '
        '0 - kmer, 1 - codon, 2 - pfam, 3 - vfam, 4 - imgHMM', required=True)
    parser.add_argument('-i', '--input',
                        help='input vector file', required=True)
    parser.add_argument('-o', '--outfile',
                        help='output file', required=True)

    args = parser.parse_args()
    file_index = args.index
    feature_code = args.feature
    file_vector = args.input
    file_output = args.outfile

    feature = Feature(file_index)
    feature_id_set = feature.get_list_to_extract(feature_code)

    vectors = VectorFile(file_vector)
    vectors.extract_feature(feature_id_set, file_output)


if __name__ == '__main__':
    main()
