#!/usr/bin/env python
from ete2 import NCBITaxa


import argparse

# Count the contigs, predicted as non-virus,  the level it got annotated
# by MEGAN. Higher score vs. lower score Are there many of them more
# ambiguous annotation

result_obj = open("all_2k_prediction.out.label", 'r')

megan_obj = open("all_2k_contigs.fa.diamond-ex-all_tax_ID.txt", 'r')

