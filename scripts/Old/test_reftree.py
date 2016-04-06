#!/usr/bin/env python
import subprocess

## Retrieve the sequence from reftree
reftreeopts = ["reftree.pl" , "--db", "genomic", "--node", "22","--attributes", "gencode"]
p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
seq = p0.stdout
print seq

