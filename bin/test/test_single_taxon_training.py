#!/usr/bin/env python
import subprocess
sttparams = ["../single_taxon_training.py", "--input", "example.mito.fasta", "--mmp", "../gm_parameters/par_11.modified","--tmpdir",".", "--taxid","12345"]
p1 = subprocess.Popen(femparams, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
metamarkout, metamarkerr= p1.communicate()
print(metamarkout)