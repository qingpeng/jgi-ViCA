#!/bin/bash 
. /global/homes/a/arrivers/dev/jgi-genelearn/bin/env.sh
module load tfmq
module load biopython
module load sqlite3

TAXON=$1
OUT=$2
CONFIG=$3

single_taxon_training.py --taxid "$TAXON" --outfile "$OUT" --config "$CONFIG"
