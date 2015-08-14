#!/usr/bin/env python


from feature_extraction_metamark import *

vect = parsemod("./test-data/")
print vect


def test_estimate_samples_lognorm_proportion():
    assert estimate_samples_lognorm(797243, 0.5, 1.333, 3000, 1140) == 69

def test_estimate_samples_fixed_proportion():
	assert estimate_samples_fixed(797243, 0.5, 5000) == 79
	
def test_estimate_samples_lognorm():
    assert estimate_samples_lognorm(797243, 200, 1.3335, 3000, 1140) == 200

def test_estimate_samples_fixed():
	assert estimate_samples_fixed(797243, 200, 5000) == 200
	
