#!/bin/bash
# copies 3 files from /global/projectb/scratch/arrivers/rna_virus/results/20150223/data/traininggenomes to
#/global/dna/projectdirs/MEP/virus/classifier/test/data/traininggenomes when executed in the first directory
# used to generate a small test set for testing create_training_data.py
FILES=`ls`
for f in $FILES
do
	cp $(find . | head -n 5) /global/dna/projectdirs/MEP/virus/classifier/test/data/traininggenomes/$f
done