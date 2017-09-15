#!/usr/bin/env python
import argparse


# 2759 - Eukaryota
# 10239 -Virus


class VectorFile:
    def __init__(self, file_vector):
        self.file_vector = file_vector

    def filter(self):
        file_vector_obj = open(self.file_vector, 'r')
        file_euk_svmlib_obj = open(self.file_vector + '.euk.svmlib', 'w')
        file_virus_svmlib_obj = open(self.file_vector + '.virus.svmlib', 'w')
        file_non_euk_svmlib_obj = open(self.file_vector + '.non_euk.svmlib',
                                       'w')
        file_non_euk_non_virus_svmlib_obj = open(
            self.file_vector + '.non_euk_non_virus.svmlib', 'w')

        file_euk_name_obj = open(self.file_vector + '.euk.name', 'w')
        file_non_euk_name_obj = open(self.file_vector + '.non_euk.name', 'w')

        file_virus_name_obj = open(self.file_vector + '.virus.name', 'w')
        file_non_euk_non_virus_name_obj = open(
            self.file_vector + '.non_euk_non_virus.name', 'w')

        n = 0
        for line in file_vector_obj:
            n += 1
            if n % 10000 == 0:
                print('...', n)
            # print line
            line = line.rstrip()
            vector_line_obj = VectorLine(line)

            tax_id = vector_line_obj.get_tax()

            if tax_id == "2759":  # euk
                file_euk_name_obj.write(vector_line_obj.get_name()+"\n")
                file_euk_svmlib_obj.write(
                    "0 "+vector_line_obj.get_svmlib()+"\n")
            else:  # non_euk
                file_non_euk_name_obj.write(vector_line_obj.get_name()+"\n")
                file_non_euk_svmlib_obj.write(
                    "1 " + vector_line_obj.get_svmlib() + "\n")

                if tax_id == "10239":  # if virus
                    file_virus_name_obj.write(vector_line_obj.get_name()+"\n")
                    file_virus_svmlib_obj.write(
                        "1 "+vector_line_obj.get_svmlib()+"\n")
                else:
                    file_non_euk_non_virus_name_obj.write(
                        vector_line_obj.get_name()+"\n")
                    file_non_euk_non_virus_svmlib_obj.write(
                        "0 "+vector_line_obj.get_svmlib()+"\n")


class VectorLine:
    def __init__(self, string_line):
        self.string_line = string_line

    def get_tax(self):
        tax_dict = {}
        #    print line
        taxonomy_list = self.string_line.split('\t')[2].split(',')[1].split(
            "=")[1].split("/")
        #    print taxonomy_list
        for taxonomy in taxonomy_list:
            groups = taxonomy.split(":")
            if groups[1] != "":
                tax_dict[groups[1]] = groups[0]
        if "0" in tax_dict:
            return tax_dict["0"]
        else:
            return "N/A"

    def get_name(self):
        name = self.string_line.split('\t')[2]
        return name

    def get_svmlib(self):
        vectors_list = self.string_line.split('\t')[3].split(" ")

        vector_strings = []
        label_id = 0
        for vector in vectors_list:
            value = float(vector.split(":")[1])
            label_id += 1
            vector_strings.append(str(label_id) + ":" + str(value))

        print_line = ' '.join(vector_strings)
        return print_line


def main():

    parser = argparse.ArgumentParser(description='A script to only keep '
                                                 'vectors about virus and '
                                                 'Eukaryota to check the '
                                                 'performance')
    parser.add_argument("vector", help='full vector file')
    args = parser.parse_args()

    vector_file_obj = VectorFile(args.vector)
    vector_file_obj.filter()

if __name__ == '__main__':
        main()
