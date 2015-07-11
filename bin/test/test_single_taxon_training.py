#!/usr/bin/env python
import subprocess
sttparams = ["single_taxon_training.py", "--input", "example.mito.fasta", "--outfile","stsoutfile.txt" , "--config", "tstt.json", "--taxid","12345"]
p1 = subprocess.Popen(sttparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
out, err= p1.communicate()
