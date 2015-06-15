#invert.py
#usage invert.py undesired_ids all_ids > output
import sys

ulines = open(sys.argv[1], 'r').read().splitlines()

with open(sys.argv[2], 'r') as all:
	for line in all:
		if not line.strip() in ulines:
			print(line.strip())