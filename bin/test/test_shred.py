#!/usr/bin/env python
import subprocess
shredparams = ["../shred.py", "--input", "example.mito.fasta", "--shred", "fixed"]
p1 = subprocess.Popen(shredparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
out, err= p1.communicate()
print(out)