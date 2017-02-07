#!/usr/bin/env python

from ete2 import NCBITaxa

file_raw_seq_input = "virus_segment_5k.fa"

file_accession2taxid_input = "prot.accession2taxid.head"

file_alignment_input = "virus_segment_5k_for_MGRAST_N.fa.daa.m8"

file_output_input = "virus_segment_5k_for_MGRAST_N.fa.daa.m8.filter"


def get_family(taxid):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    ranks = ncbi.get_rank(lineage)
    for rank in ranks.keys():
        if ranks[rank] == 'family':
            return rank

    return "None"


def test_same_family(taxid1, taxid2):
    family1 = get_family(taxid1)
    family2 = get_family(taxid2)
    # if any one does not have family information, keep it
    if family1 == "None" or family2 == "None":
        return False

    if family1 == family2:
        return True
    else:
        return False


def main(file_raw_seq, file_accession2taxid, file_alignment, file_output):

    file_accession2taxid_obj = open(file_accession2taxid, 'r')

    dict_accession2taxid = {}
    file_accession2taxid_obj.readline()
    for line in file_accession2taxid_obj:
        line = line.rstrip()
        fields = line.split()
        dict_accession2taxid[fields[0]] = int(fields[2])

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
    file_output_obj = open(file_output, 'w')

    for line in file_alignment_obj:
        line = line.rstrip()
        fields = line.split()
        query = int(fields[0])
        target = fields[1].split(".")[0]
        taxid_query = dict_seqid_taxid[query]
        taxid_target = dict_accession2taxid[target]

        if not test_same_family(taxid_query, taxid_target):
            file_output_obj.write(line+'\n')


main(file_raw_seq_input, file_accession2taxid_input, file_alignment_input, file_output_input)





