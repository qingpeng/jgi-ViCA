
# A script to take a directory of metamark models and create a profile vector
#metamarkparser directory class
import os
import sys
import re
from itertools import chain

directory = sys.argv[1]
cls = sys.argv[2]

def flatten(listOfLists):
    """Flatten one level of nesting"""
    return chain.from_iterable(listOfLists)

for file in os.listdir(directory):
	if file.endswith("hmm.mod"):
		with open(directory +'/'+file) as f:
			rawvec = []
			for line in f.readlines():
				rawvec.append(line.split())
			#print(compvec)
			start1 = rawvec.index(['COD1']) +1
			stop1 = rawvec.index(['COD2']) -1
			COD1 = (list(flatten(rawvec[start1: stop1 ])))
			start2 = 1 + rawvec.index(['NONC'])
			NONC = (list(flatten(rawvec[start2: start2 +64 ])))
			print(cls+'\t'+'\t'.join(COD1+NONC))
        		
        		
