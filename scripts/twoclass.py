#!/usr/bin/env python
#twoclass.py
import sys

with open(sys.argv[1], 'r') as file:
	for line in file:
		ll = line.split('\t')
		if ll[0] in ["ssRNAphage","ssRNAvirus", "ssDNAphage","ssDNAvirus","Retroviruses", "dsRNAphage", "dsRNAvirus", "dsDNAvirus", "dsDNAphage"]:
			ll[0] = "viral"
		else:
			ll[0] = "nonviral"
		print('\t'.join(ll))
file.close()