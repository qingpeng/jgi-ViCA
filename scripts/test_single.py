#!/usr/bin/env python
import subprocess


shredopts = ["shred.py",  "--shred", "fixed", "--samples", "200", \
    "--length", "5000", "--output", "shred2.fa","--input","138_genome.fa"]
    
p1 = subprocess.Popen(shredopts, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
