#!/usr/bin/env python

import argparse


class SVMLIB:
    def __init__(self, file_name):
        self.file_name = file_name

    def pick(self, output_name):
        file_in_obj = open(self.file_name, 'r')
        file_out_obj = open(output_name, 'w')
        for line in file_in_obj:
            line = line.rstrip()
            fields = line.split()
            if len(fields) >= 393:
                new_line = ' '.join(fields[0:393])
            else:
                new_line = line
            file_out_obj.write(new_line+'\n')


def main():
    parser = argparse.ArgumentParser(
        description='A script to use remove pfam features from svmlib file')

    parser.add_argument('libsvm', help='libsvm file of training data')
    parser.add_argument('output',  help='output file')

    args = parser.parse_args()
    svmlib = SVMLIB(args.libsvm)
    svmlib.pick(args.output)

if __name__ == '__main__':
    main()
