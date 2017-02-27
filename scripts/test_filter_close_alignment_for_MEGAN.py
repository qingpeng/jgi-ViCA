from .filter_close_alignment_for_MEGAN \
    import AlignmentFile, AccessionTaxidFile, SequenceFile, TaxID


def test_taxid_get_rank():
    tax_id_1 = TaxID('10090')
    assert tax_id_1.get_rank('family') == 10066
    assert tax_id_1.get_rank('famiy') == "N/A"
    tax_id_1 = TaxID('10091')
    assert tax_id_1.get_rank('genus') == 10088


def test_accession_taxid_file():
    accession_taxid_file = AccessionTaxidFile(
        '../tests/test-data/prot.accession2taxid.head')
    target_set1 = {'WP_072633130', 'P06912', 'P06913'}
    dict_accession2taxid = accession_taxid_file.get_dict_accession2taxid(
        target_set1)
    assert dict_accession2taxid == {'WP_072633130': 10090, 'P06912': 9986}


def test_sequence_file_get_dict_seqid_taxid():
    sequence_file = SequenceFile('../tests/test-data/test_segment_5k.fa.3')
    dict_seqid_taxid = sequence_file.get_dict_seqid_taxid()
    assert dict_seqid_taxid == {1: 101570, 2: 101570, 3: 1029988}


def test_alignment_file_get_target_set():
    alignment_file = AlignmentFile(
        '../tests/test-data/test.m8', '../tests/test-data/test.m8.filter',
        '../tests/test-data/test_segment_5k.fa.3',
        '../tests/test-data/prot.accession2taxid.head', 3)
    target_set = alignment_file.get_target_set(alignment_file.file_alignment)
    assert target_set == {'WP_049883809', 'WP_060456234', 'WP_015976792',
                          'KWR46313', 'WP_024156024', 'WP_047064325',
                          'ETE43525', 'WP_070818888', 'WP_063889886',
                          'WP_072633130'}


def test_alignment_file_test_same_rank():
    assert not AlignmentFile.test_same_rank(10291, 10591, 'genus')
    assert not AlignmentFile.test_same_rank(10091, 10591, 'order')
    assert AlignmentFile.test_same_rank(10090, 10091, 'order')
    assert not AlignmentFile.test_same_rank(10090, 10291, 'family')


def test_alignment_file_get_query_target():
    line = '1	WP_070818888.1	76.8	924	214	0	3	' \
           '2774	201	1124	0.0e+00	1239.2'
    query, target = AlignmentFile.get_query_target(line)
    assert query == 1
    assert target == 'WP_070818888'


def test_alignment_file_check_line():
    alignment_file = AlignmentFile(
        '../tests/test-data/test.m8', '../tests/test-data/test.m8.filter',
        '../tests/test-data/test_segment_5k.fa.3',
        '../tests/test-data/prot.accession2taxid.head', 3)
    line = '1	WP_070818888.1	76.8	924	214	0	3	' \
           '2774	201	1124	0.0e+00	1239.2'
    alignment_file.dict_accession2taxid = {
        'WP_072633130': 10090, 'P06912': 9986}
    alignment_file.dict_seqid_taxid = {1: 10090, 2: 101570, 3: 1029988}
    assert alignment_file.check_line(line, 'family')  # access does not exit

    alignment_file.dict_seqid_taxid = {1: 10090, 2: 101570, 3: 1029988}
    alignment_file.dict_accession2taxid = {
        'WP_070818888': 10591, 'P06912': 9986}
    assert alignment_file.check_line(line, 'order')

    alignment_file.dict_seqid_taxid = {1: 10090, 2: 101570, 3: 1029988}
    alignment_file.dict_accession2taxid = {
        'WP_070818888': 10091, 'P06912': 9986}
    assert not alignment_file.check_line(line, 'order')


def test_alignment_file_filtering():
    alignment_file = AlignmentFile(
        '../tests/test-data/test.m8', '../tests/test-data/test.m8.filter',
        '../tests/test-data/test_segment_5k.fa.3',
        '../tests/test-data/prot.accession2taxid.head', 3)

    alignment_file.dict_accession2taxid = {
        'WP_015976792': 10091, 'WP_070818888': 10091, 'WP_072633130': 10091,
        'ETE43525': 10090, 'WP_049883809': 10090, 'WP_024156024': 10090,
        'KWR46313': 10291, 'WP_047064325': 10591, 'WP_063889886': 10791,
        'WP_060456234': 10796}
    alignment_file.dict_seqid_taxid = {1: 10090, 2: 10111, 3: 10090}
    assert alignment_file.filtering() == 1






#
# def test_taxid():
#     tax_id_1 = TaxID('')
#     print dict_seqid_taxid
#     assert dict_seqid_taxid == {}

