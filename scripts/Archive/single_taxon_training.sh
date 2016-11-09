#!/bin/bash 
. /global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/env.sh

TAXON=$1
OUT=$2
CONFIG=$3

single_taxon_training.py --taxid "$TAXON" --outfile "$OUT" --config "$CONFIG"
