#!/usr/bin/env python

# 02/07 can do the filtering for all 4 levels altogether now


from ete2 import NCBITaxa
import sys

import argparse


def get_rank(taxid, rank_level):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    for rank in ranks.keys():
        if ranks[rank] == rank_level:
            return rank

    return "None"


def test_same_rank(taxid1, taxid2, rank_level):
    family1 = get_rank(taxid1, rank_level)
    family2 = get_rank(taxid2, rank_level)
    # print family1, family2
    # if any one does not have family information, keep it
    if family1 == "None" or family2 == "None":
        # print "false"
        return False

    if family1 == family2:
        # print "true"
        return True
    else:
        # print "false"
        return False


def filtering(file_raw_seq, file_accession2taxid, file_alignment,
              file_output, rank_level, top_number):
    file_alignment_obj = open(file_alignment, 'r')

    target_set = set()
    for line in file_alignment_obj:
        line = line.rstrip()
        fields = line.split()
        target = fields[1].split(".")[0]
        target_set.add(target)

    print "file_alignment done!"

    file_accession2taxid_obj = open(file_accession2taxid, 'r')

    dict_accession2taxid = {}
    file_accession2taxid_obj.readline()
    count = 0
    block = 0
    num_hit = {}
    for line in file_accession2taxid_obj:
        # print line
        count = count + 1
        if count == 1000000:
            block += 1
            count = 0
            print "1111111111111111111111111", block, "\n"
        line = line.rstrip()
        fields = line.split()
        if fields[0] in target_set:
            # print "hit!"
            dict_accession2taxid[fields[0]] = int(fields[2])

    print "file_accession2taxid done!"

    file_raw_seq_obj = open(file_raw_seq, 'r')
    seqid = 0
    dict_seqid_taxid = {}
    for line in file_raw_seq_obj:
        line = line.rstrip()
        if line[0] == ">":
            taxid = line.split("taxid=")[1].split(",")[0]
            seqid += 1
            dict_seqid_taxid[seqid] = int(taxid)

    file_alignment_obj = open(file_alignment, 'r')

    if rank_level == "all":
        file_output_family_obj = open(file_output+'.family', 'w')
        file_output_order_obj = open(file_output+'.order', 'w')
        file_output_genus_obj = open(file_output+'.genus', 'w')
        file_output_species_obj = open(file_output+'.species', 'w')
    else:
        file_output_obj = open(file_output, 'w')

    block = 0
    count = 0
    for line in file_alignment_obj:
        line = line.rstrip()
        fields = line.split()
        query = int(fields[0])
        target = fields[1].split(".")[0]

        count = count + 1
        if count == 1000000:
            block += 1
            count = 0
            print "22222222", block, "\n"

        if num_hit.setdefault(query, 0) == top_number:
            continue

        try:  # make sure they belong to different family explicitly,
            # if no family info, remove the hit

            taxid_query = dict_seqid_taxid[query]
            taxid_target = dict_accession2taxid[target]
            # print "test", taxid_query, taxid_target

            if rank_level == "all":
                if not test_same_rank(taxid_query, taxid_target, "family"):
                    # print "beforePprint"
                    num_hit[query] += 1
                    # print "print"
                    file_output_family_obj.write(line + '\n')

                if not test_same_rank(taxid_query, taxid_target, "order"):
                    # print "beforePprint"
                    num_hit[query] += 1
                    # print "print"
                    file_output_order_obj.write(line + '\n')

                if not test_same_rank(taxid_query, taxid_target, "genus"):
                    # print "beforePprint"
                    num_hit[query] += 1
                    # print "print"
                    file_output_genus_obj.write(line + '\n')

                if not test_same_rank(taxid_query, taxid_target, "species"):
                    # print "beforePprint"
                    num_hit[query] += 1
                    # print "print"
                    file_output_species_obj.write(line + '\n')

            else:

                if not test_same_rank(taxid_query, taxid_target, rank_level):
                    # print "beforePprint"
                    num_hit[query] += 1
                    # print "print"
                    file_output_obj.write(line + '\n')
        except KeyError:
            continue

    if rank_level == "all":
        file_output_family_obj.close()
        file_output_order_obj.close()
        file_output_genus_obj.close()
        file_output_species_obj.close()
    else:
        file_output_obj.close()

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
    parser.add_argument('-r', '--rank',
                        help="rank level all,family,order", required=True)
    parser.add_argument('-t', '--top',
                        help="number of top hits to keep, 100", required=True)
    # virus_segment_5k_for_MGRAST_N.fa.daa.m8.filter
    parser.add_argument('-o', '--output',
                        help="filtered alignment output", required=True)

    args = parser.parse_args()

    file_raw_seq_input = args.sequence  # "virus_segment_5k.fa"

    file_accession2taxid_input = args.accession  # "prot.accession2taxid"

    file_alignment_input = args.input
    # "virus_segment_5k_for_MGRAST_N.fa.daa.m8"

    file_output_input = args.output
    # "virus_segment_5k_for_MGRAST_N.fa.daa.m8.filter"

    rank_level = args.rank  # "family" or "order" or "all"

    if rank_level not in ['all', 'family', 'order', 'genus', 'species']:
        print "Error: rank_level wrong input, must be all, " \
              "family, order, genus, species"
        exit()

    top_number = int(args.top)

    filtering(file_raw_seq_input,
              file_accession2taxid_input, file_alignment_input,
              file_output_input, rank_level, top_number)


main()
