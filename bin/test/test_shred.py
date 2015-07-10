#!/usr/bin/env python
import subprocess

print("Testing shred.py with fixed length")
shredparams = ["../shred.py", "--input", "example.mito.fasta", "--shred", "fixed"]
p1 = subprocess.Popen(shredparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
out, err= p1.communicate()
print(out)

print("Testing shred.py with length set by a gamma distribution")
shredparams2 = ["../shred.py", "--input", "example.mito.fasta", "--shred", "gamma", "--alpha", "5","--scale", "45", "--loc", "100" ]
p2 = subprocess.Popen(shredparams2, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
out2, err2= p2.communicate()
print(out2)