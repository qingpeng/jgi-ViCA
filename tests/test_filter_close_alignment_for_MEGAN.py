from scripts.filter_close_alignment_for_MEGAN \
    import AlignmentFile, AccessionTaxidFile, SequenceFile, TaxID


def test_taxid_get_rank():
    tax_id_1 = TaxID('10090')
    assert tax_id_1.get_rank('family') == 10066
    assert tax_id_1.get_rank('famiy') == "N/A"
    tax_id_1 = TaxID('10091')
    assert tax_id_1.get_rank('genus') == 10088


def test_accession_taxid_file():
    accession_taxid_file = AccessionTaxidFile(
        './test-data/prot.accession2taxid.head')
    target_set1 = {'WP_015976792', 'ETE43524', 'WP_075180623'}
    dict_accession2taxid = accession_taxid_file.get_dict_accession2taxid(
        target_set1)
    assert dict_accession2taxid == {'WP_015976792': 9606,
                                    'WP_075180623': 1925589}


def test_sequence_file_get_dict_seqid_taxid():
    sequence_file = SequenceFile('../tests/test-data/test_segment_5k.fa.3')
    dict_seqid_taxid = sequence_file.get_dict_seqid_taxid()
    assert dict_seqid_taxid == {1: 101570, 2: 10090, 3: 10091}


def test_alignment_file_get_target_set():
    alignment_file = AlignmentFile(
        './test-data/test.m8', './test-data/test.m8.filter',
        './test-data/test_segment_5k.fa.3',
        './test-data/prot.accession2taxid.head', 3)
    target_set = alignment_file.get_target_set(alignment_file.file_alignment)
    print target_set
    assert target_set == {'WP_049883809', 'WP_015976792', 'WP_075180623',
                          'KWR46313', 'WP_024156024', 'WP_024156026',
                          'ETE43525', 'WP_070818888', 'WP_072633130'}


def test_alignment_file_test_same_rank_keep_ambiguous():
    # 10291 - no genus 10591 - "333750" - return False, not same rank
    assert not AlignmentFile.test_same_rank_keep_ambiguous(10291, 10591,
                                                           'genus')
    # 10091 - 9989, 10591 - no order - not same rank
    assert not AlignmentFile.test_same_rank_keep_ambiguous(10091, 10591,
                                                           'order')
    # 10090 - 9989, 10091 - 9989, return True
    assert AlignmentFile.test_same_rank_keep_ambiguous(10090, 10091, 'order')

    # 10090 -10066, 10291 - 10240, not the same, return False
    assert not AlignmentFile.test_same_rank_keep_ambiguous(10090, 10291,
                                                           'family')
    # 101570 - 10699 family, 1925589 - 1903409, return False
    assert not AlignmentFile.test_same_rank_keep_ambiguous(101570, 1925589,
                                                           'family')


def test_alignment_file_test_same_rank_discard_ambiguous():
    # 10291 - no genus 10591 - "333750" - return True, treat as same rank
    assert AlignmentFile.test_same_rank_discard_ambiguous(10291, 10591,
                                                          'genus')
    # 10091 - 9989, 10591 - no order - return True, treat as same rank
    assert AlignmentFile.test_same_rank_discard_ambiguous(10091, 10591,
                                                          'order')
    # 10090 - 9989, 10091 - 9989, return True
    assert AlignmentFile.test_same_rank_discard_ambiguous(10090, 10091,
                                                          'order')
    # 10090 -10066, 10291 - 10240, not the same, return False
    assert not AlignmentFile.test_same_rank_discard_ambiguous(10090, 10291,
                                                              'family')
    # 101570 - 10699 family, 1925589 - 1903409, return False
    assert not AlignmentFile.test_same_rank_keep_ambiguous(101570, 1925589,
                                                           'family')


def test_alignment_file_get_query_target():
    line = '1	WP_070818888.1	76.8	924	214	0	3	' \
           '2774	201	1124	0.0e+00	1239.2'
    query, target = AlignmentFile.get_query_target(line)
    assert query == 1
    assert target == 'WP_070818888'


def test_alignment_file_check_line():
    alignment_file = AlignmentFile(
        './test-data/test.m8', './test-data/test.m8.filter',
        './test-data/test_segment_5k.fa.3',
        './test-data/prot.accession2taxid.head', 3)
    line = '1	WP_070818888.1	76.8	924	214	0	3	' \
           '2774	201	1124	0.0e+00	1239.2'

    # target accession does not exist
    alignment_file.dict_accession2taxid = {
        'WP_072633130': 10090, 'P06912': 9986}
    alignment_file.dict_seqid_taxid = {1: 10090, 2: 101570, 3: 1029988}
    assert not alignment_file.check_line(
        line, 'family', 'discard')

    alignment_file.dict_seqid_taxid = {1: 10090, 2: 101570, 3: 1029988}
    alignment_file.dict_accession2taxid = {
        'WP_070818888': 10591, 'P06912': 9986}
    assert not alignment_file.check_line(
        line, 'order', 'keep')
    assert alignment_file.check_line(
        line, 'order', 'discard')

    line = '1	WP_075180623.1	41.6	735	415	7	588	2774	320	1046	' \
           '8.4e-134	488.0'
    alignment_file.dict_seqid_taxid = {1: 101570, 2: 10090, 3: 1029988}
    alignment_file.dict_accession2taxid = {
        'WP_070818888': 10091, 'WP_075180623': 1925589}
    assert not alignment_file.check_line(
        line, 'family', 'discard')
    assert not alignment_file.check_line(
        line, 'family', 'keep')


def test_alignment_file_filtering_keep_ambiguous():
    alignment_file = AlignmentFile(
        './test-data/test.m8', './test-data/test.m8.filter.keep',
        './test-data/test_segment_5k.fa.3',
        './test-data/prot.accession2taxid.head', 3)
    # m8, 1: 101570, 2: 10090, 3: 10091
    assert alignment_file.filtering_with_option('keep') == 1


def test_alignment_file_filtering_discard_ambiguous():
    alignment_file = AlignmentFile(
        './test-data/test.m8',
        './test-data/test.m8.filter.discard',
        './test-data/test_segment_5k.fa.3',
        './test-data/prot.accession2taxid.head', 3)

    assert alignment_file.filtering_with_option('discard') == 1

