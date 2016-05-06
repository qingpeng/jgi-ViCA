#!/usr/bin/env python
#twoclass.py
import sys

with open(sys.argv[1], 'r') as file:
	for line in file:
		ll = line.split('\t')
#		print ll
		if ll[0] in ["ssRNAPhage","ssRNAVirus", "dsRNAPhage","dsRNAVirus","dsDNAPhage", "dsDNAVirus", "ssDNAPhage", "ssDNAVirus", "Retroviruses"]:
			ll[0] = "viral"
#			print "viral"
		else:
			ll[0] = "nonviral"
#			print "nonviral"
		print('\t'.join(ll))
file.close()