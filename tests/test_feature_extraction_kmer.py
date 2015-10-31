#!/usr/bin/env python


from feature_extraction_kmer import *

def test_parsemod():
    vect = parsemod("./test-data/")
    assert len(vect) == 256
    assert vect[191] == '0.015209'
    assert vect[193] == '0.016226'
    assert vect[-1] == '0.016415'

