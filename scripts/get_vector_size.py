#!/usr/bin/env python
import argparse


class VectorFile:
    def __init__(self, vector_file):
        self.vector_file = vector_file

    def get_size(self):
        file_vector_file_obj = open(self.vector_file, 'r')
        file_output_obj = open(self.vector_file+'.size', 'w')
        max_id = 0
        count = 0
        for line in file_vector_file_obj:
            line = line.rstrip()
            fields = line.split()
            feature_id = int(fields[-1].split(":")[0])
            if feature_id > max_id:
                max_id = feature_id
            count += 1
        line_output1 = 'number of rows: '+str(count)
        line_output2 = 'size of vector: '+str(max_id)
        file_output_obj.write(line_output1+'\n')
        file_output_obj.write(line_output2+'\n')


def main():
    parser = argparse.ArgumentParser(
        description='A script to count the size of vector file')
    parser.add_argument('-i', '--input',
                        help='input vector file', required=True)

    args = parser.parse_args()
    file_vector = args.input

    vector_file = VectorFile(file_vector)
    vector_file.get_size()

if __name__ == '__main__':
    main()
