#!/bin/bash 
. /global/homes/q/qpzhang/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/env.sh

#module load tfmq
module load biopython
module load sqlite3
module load hmmer/3.1b2

TAXON=$1
OUT=$2
CONFIG=$3

single_taxon_training.py --taxid "$TAXON" --outfile "$OUT" --config "$CONFIG"
