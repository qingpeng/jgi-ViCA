#!/usr/bin/env python

# 02/07 can do the filtering for all 4 levels altogether now


from ete2 import NCBITaxa


import argparse


class TaxID:
    def __init__(self, tax_id):
        self.tax_id = tax_id

    def get_rank(self, rank_level):
        ncbi = NCBITaxa()
        lineage = ncbi.get_lineage(self.tax_id)
        ranks = ncbi.get_rank(lineage)
        for rank in ranks.keys():
            if ranks[rank] == rank_level:
                return rank

        return "N/A"


class AccessionTaxidFile:
    def __init__(self, file_name):
        self.file_name = file_name
        self.dict_accession2taxid = {}

    def get_dict_accession2taxid(self, target_set):
        file_accession2taxid_obj = open(self.file_name, 'r')
        file_accession2taxid_obj.readline()
        for line in file_accession2taxid_obj:
            # print line
            line = line.rstrip()
            fields = line.split()
            if fields[0] in target_set:
                # print "hit!"
                self.dict_accession2taxid[fields[0]] = int(fields[2])
        file_accession2taxid_obj.close()
        print "AccessionTax file processing done!"
        return self.dict_accession2taxid


class SequenceFile:
    def __init__(self, file_name):
        self.file_name = file_name
        self.dict_seqid_taxid = {}

    def get_dict_seqid_taxid(self):

        file_raw_seq_obj = open(self.file_name, 'r')
        seqid = 0
        for line in file_raw_seq_obj:
            line = line.rstrip()
            if line[0] == ">":
                taxid = line.split("taxid=")[1].split(",")[0]
                seqid += 1
                self.dict_seqid_taxid[seqid] = int(taxid)
        file_raw_seq_obj.close()
        print "Reading fasta file done!"
        return self.dict_seqid_taxid


class AlignmentFile:
    def __init__(self, file_alignment, file_output, file_seq, file_accession,
                 top_number):
        self.file_alignment = file_alignment
        self.file_output = file_output
        self.file_seq = file_seq
        self.file_accession = file_accession
        self.top_number = top_number

        self.dict_seqid_taxid = SequenceFile(
            self.file_seq).get_dict_seqid_taxid()
        self.target_set = self.get_target_set(self.file_alignment)
        self.dict_accession2taxid = AccessionTaxidFile(
            self.file_accession).get_dict_accession2taxid(self.target_set)

    @staticmethod
    def get_target_set(file_alignment):
        file_alignment_obj = open(file_alignment, 'r')
        target_set = set()
        for line in file_alignment_obj:
            line = line.rstrip()
            fields = line.split()
            target = fields[1].split(".")[0]  # remove version
            target_set.add(target)
        file_alignment_obj.close()
        print "Get target set done!"
        return target_set

    @staticmethod
    def test_same_rank_keep_ambiguous(
            tax_id_1, tax_id_2, rank_level):
        tax_id_1_rank = TaxID(tax_id_1).get_rank(rank_level)
        tax_id_2_rank = TaxID(tax_id_2).get_rank(rank_level)
        # if any one does not have family information, keep it, not the same!
        if tax_id_1_rank == "N/A" or tax_id_2_rank == "N/A":
            return False

        if tax_id_1_rank == tax_id_2_rank:
            return True
        else:  # only if the tax_id on that level does not equal explicitly
            return False

    @staticmethod
    def test_same_rank_discard_ambiguous(
            tax_id_1, tax_id_2, rank_level):
        tax_id_1_rank = TaxID(tax_id_1).get_rank(rank_level)
        tax_id_2_rank = TaxID(tax_id_2).get_rank(rank_level)

        # if any one does not have family information, discard it,
        # can't make
        # sure it is in the same family or not, probably it is!
        if tax_id_1_rank == "N/A" or tax_id_2_rank == "N/A":
            return True

        if tax_id_1_rank == tax_id_2_rank:
            return True
        else:  # only if the tax_id on that level does not equal explicitly
            return False

    @staticmethod
    def get_query_target(line):
        line = line.rstrip()
        fields = line.split()
        query = int(fields[0])
        target = fields[1].split(".")[0]
        return query, target

    def check_line(self, line, tax_level, ambiguous_option):
        """
        check if the line of hit has query and target in the same tax group
        :param line:
        :param tax_level:
        :param ambiguous_option: "keep" or "discard"
        :return: True if query and target in the same tax group
                 False if not
        """
        query, target = self.get_query_target(line)
        taxid_query = self.dict_seqid_taxid[query]
        try:
            taxid_target = self.dict_accession2taxid[target]
        # print "test", taxid_query, taxid_target
        except KeyError:
            print "Can't find taxid for target:", target
            return False  # if can't find taxid for target, don't keep
            #  this record.
        # print taxid_query, taxid_target, tax_level
        # print self.test_same_rank(taxid_query, taxid_target, tax_level)
        if ambiguous_option == "keep":
            return AlignmentFile.test_same_rank_keep_ambiguous(
                    taxid_query, taxid_target, tax_level)
        else:
            return AlignmentFile.test_same_rank_discard_ambiguous(
                    taxid_query, taxid_target, tax_level)

    def filtering_with_option(self, ambiguous_option):

        file_alignment_obj = open(self.file_alignment, 'r')

        file_output_family_obj = open(self.file_output+'.family', 'w')
        file_output_order_obj = open(self.file_output+'.order', 'w')
        file_output_genus_obj = open(self.file_output+'.genus', 'w')

        num_hit_family = {}
        num_hit_order = {}
        num_hit_genus = {}

        for line in file_alignment_obj:
            line = line.rstrip()
            query, target = self.get_query_target(line)

            if query not in num_hit_genus:
                num_hit_genus[query] = 0
            if query not in num_hit_family:
                num_hit_family[query] = 0
            if query not in num_hit_order:
                num_hit_order[query] = 0

            if (num_hit_genus[query] == self.top_number
                and num_hit_family[query] == self.top_number
                    and num_hit_order[query] == self.top_number):

                continue

            # make sure they belong to different family explicitly,
            # if no family info, remove the hit

            if num_hit_family[query] != self.top_number:
                if not self.check_line(line, "family", ambiguous_option):
                    # print "beforePprint"
                    num_hit_family[query] += 1
                    # print "print"
                    file_output_family_obj.write(line + '\n')

            if num_hit_order[query] != self.top_number:
                if not self.check_line(line, "order", ambiguous_option):
                    # print "beforePprint"
                    num_hit_order[query] += 1
                    # print "print"
                    file_output_order_obj.write(line + '\n')

            if num_hit_genus[query] != self.top_number:
                if not self.check_line(line, "genus", ambiguous_option):
                    # print "beforePprint"
                    num_hit_genus[query] += 1
                    # print "print"
                    file_output_genus_obj.write(line + '\n')

        file_output_family_obj.close()
        file_output_order_obj.close()
        file_output_genus_obj.close()
        return 1


def main():
    parser = argparse.ArgumentParser(
        description='A script to filter out the alignments hits close to the '
                    'query in specific taxonomy rank level, or all taxonomy '
                    'rank levels, - family, order, genus, species')
    parser.add_argument('-s', '--sequence',
                        help="raw sequence,virus_segment_5k.fa", required=True)
    parser.add_argument('-a', '--accession',
                        help="prot.accession2taxid", required=True)
    parser.add_argument('-i', '--input',
                        help="alignment file to filter, .m8 format",
                        required=True)
    parser.add_argument('-t', '--top',
                        help="number of top hits to keep, 100", required=True)
    # virus_segment_5k_for_MGRAST_N.fa.daa.m8.filter
    parser.add_argument('-o', '--output',
                        help="filtered alignment output", required=True)
    parser.add_argument('--filter_option',
                        help="option for filtering, "
                             "discard - discard alignment if query/target does"
                             "not have tax_id on specific level, "
                             "keep - keep those alignments", required=True)

    args = parser.parse_args()

    file_raw_seq_input = args.sequence  # "virus_segment_5k.fa"
    file_accession2taxid_input = args.accession  # "prot.accession2taxid"
    file_alignment_input = args.input
    # "virus_segment_5k_for_MGRAST_N.fa.daa.m8"
    file_output_input = args.output
    # "virus_segment_5k_for_MGRAST_N.fa.daa.m8.filter"
    top_number = int(args.top)
    filter_option = str(args.filter_option)

    alignment_file_obj = AlignmentFile(
        file_alignment_input, file_output_input, file_raw_seq_input,
        file_accession2taxid_input, top_number)
    alignment_file_obj.filtering_with_option(filter_option)


if __name__ == '__main__':
    main()
