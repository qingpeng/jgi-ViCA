#!/usr/bin/env python
import argparse
import math


def fun_log(value):
    if value == 0.0:
        value = "1e-200"
    return math.log(float(value))


def get_feature_list(file_feature):
    """ Get feature list from the file with feature list"""
    file_feature_obj = open(file_feature, 'r')
    feature_list = []
    for line in file_feature_obj:
        line = line.rstrip()
        fields = line.split()
        feature_list.append(fields[1])
    return feature_list


def to_svmlib(vectors_line, feature_list):
    """ convert vector line into svmlib format
    keep adding new label into feature_list
    """

    vectors_list = vectors_line.split(" ")

    vector_strings = []

    for vector in vectors_list:
        label = vector.split(":")[0]
        pre = label.split("_")[0]
        value = float(vector.split(":")[1])
        if pre == "2" or pre == "3" or pre == "4":
            # if pfam or vfam or IMG vfam ,get log
            value = fun_log(value)

        if label in feature_list:
            label_id = feature_list.index(label) + 1
            vector_strings.append(str(label_id) + ":" + str(value))

    dict_vectors = {}
    for vector in vector_strings:
        label = vector.split(":")
        dict_vectors[int(label[0])] = label[1]
    # sort the feature id
    sorted_key = sorted(dict_vectors.keys())
    print sorted_key
    out = ''
    for key in sorted_key:
        out = out + ' ' + str(key) + ':' + dict_vectors[key]

    return out


def convert_file(file_vector, file_output, file_feature):
    """ main function """
    feature_list = get_feature_list(file_feature)
    # get feature list used in training, only keep the feature in this list
    file_vector_obj = open(file_vector, 'r')
    file_output_obj = open(file_output, 'w')
    for line in file_vector_obj:
        line = line.rstrip()
        fields = line.split('\t')
        vectors_line = fields[3]
        svmlib_line = to_svmlib(vectors_line, feature_list)
        print_line = "0" + svmlib_line  # always label as 0 here
        file_output_obj.write(print_line+'\n')


def main():
    parser = argparse.ArgumentParser(
        description='A script to convert vector file to libsvm for prediction')

    parser.add_argument('vector', help='full vector file')
    parser.add_argument('feature_file',
                        help='feature listfile, from training set, with model')
    parser.add_argument('outfile', help='output file')

    args = parser.parse_args()
    convert_file(args.vector, args.outfile, args.feature_file)

if __name__ == '__main__':
    main()
