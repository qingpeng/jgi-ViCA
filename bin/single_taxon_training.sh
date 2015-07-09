#!/bin/bash
. env.sh
#. /global/homes/e/ekirton/dev/RefTree/env.sh

module load biopython

TAXON=$1
OUT=$2
python single_taxon_training.py --taxid '$TAXON' --outfile '$OUT' 
