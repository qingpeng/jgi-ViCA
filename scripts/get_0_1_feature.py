#!/usr/bin/env python
import argparse



class VectorLine:
    def __init__(self, string_line):
        self.string_line = string_line

    def get_svmlib(self):
        vectors_list = self.string_line.split(' ')
        return ' '.join(vectors_list[0:393])



def main():

    parser = argparse.ArgumentParser(description='A script to only keep '
                                                 'vectors about virus and '
                                                 'Eukaryota to check the '
                                                 'performance')
    parser.add_argument("vector", help='full vector file')
    args = parser.parse_args()

    file_in_obj = open(args.vector, 'r')
    file_out_obj = open(args.vector+'.0_1', 'w')

    for line in file_in_obj:
        line.rstrip()
        vectorline = VectorLine(line)
        outline = vectorline.get_svmlib()
        file_out_obj.write(outline+'\n')



if __name__ == '__main__':
        main()