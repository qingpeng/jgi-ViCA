#!/usr/bin/env python
import argparse


class Hitsfile:
    def __init__(self, file_in):
        self.file_in = file_in

    def get_refseq(self, output_file):
        file_inputfile_obj = open(self.file_in, 'r')
        file_output_obj = open(output_file, 'w')

        for line in file_inputfile_obj:
            line = line.rstrip()
            fields = line.split()
            prefix = fields[1][:3]
            if (prefix == "NP" or prefix == 'AP_' or prefix == 'XP_' or
                    prefix == 'YP_' or prefix == 'WP_'):
                file_output_obj.write(line+'\n')
        file_output_obj.close()


def main():
    parser = argparse.ArgumentParser(
        description='A script to extract refseq related hits from m8 file')
    parser.add_argument('-i', '--input',
                        help='input m8 file', required=True)
    parser.add_argument('-o', '--output',
                        help='output m8 file', required=True)

    args = parser.parse_args()
    file_in = args.input
    file_out = args.output

    hits_file = Hitsfile(file_in)
    hits_file.get_refseq(file_out)


if __name__ == '__main__':
    main()
