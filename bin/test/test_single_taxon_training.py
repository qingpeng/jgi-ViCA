#!/usr/bin/env python
import subprocess
sttparams = ["single_taxon_training.py", "--outfile","stsoutfile.txt" , "--config", "tstt.json", "--taxid","246200"]
p1 = subprocess.Popen(sttparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err= p1.communicate()
print(out)
print(err)
