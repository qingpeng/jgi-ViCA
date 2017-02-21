#!/usr/bin/env python
import argparse
import random
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

        try:
            label_id = feature_list.index(label)
        except ValueError:
            feature_list.append(label)
            label_id = feature_list.index(label)

        vector_strings.append(str(label_id) + ":" + str(value))

    print_line = ' '.join(vector_strings)
    return print_line, feature_list


def write_new_feature_file(feature_list, new_file_feature):
    """
    :param feature_list:
    :param new_file_feature:
    :return: just write new feature list into the file
    """
    file_new_feature_obj = open(new_file_feature, 'w')
    count = 0
    for feature in feature_list:
        line = str(count) + ' ' + feature
        file_new_feature_obj.write(line + '\n')
        count += 1
    file_new_feature_obj.close()


def convert_file(file_vector,file_output,file_feature,target_label):
    """ main function """
    feature_list = get_feature_list(file_feature)
    file_vector_obj = open(file_vector, 'r')
    file_output_obj = open(file_output, 'w')
    for line in file_vector_obj:
        line = line.rstrip()
        fields = line.split('\t')
        vectors_line = fields[3]
        svmlib_line, feature_list = to_svmlib(vectors_line, feature_list)
        print_line = target_label + ' ' + svmlib_line
        file_output_obj.write(print_line+'\n')
    new_file_feature = file_feature + '.new_list'
    write_new_feature_file(feature_list, new_file_feature)


def main():
    parser = argparse.ArgumentParser(
        description='A script to convert vector file to svmlib')
    parser.add_argument('-v', '--vector', help='full vector file',
                        required=True)
    parser.add_argument('-f', '--feature_file',
                        help='feature list file, or None',
                        required=True)
    parser.add_argument('-o', '--outfile', help='output file', required=True)

    parser.add_argument('-t', '--target',
                        help='target label to assign', required=True)

    args = parser.parse_args()


    convert_file(args.vector, args.outfile, args.feature_file, args.target)


if __name__ == '__main__':
    main()


