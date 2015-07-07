#!/bin/sh

. env.sh
module load biopython
#Designed to be submitted to 8 core commodity nodes with two threads per worker
exec tfmq-worker -z -q singletaxontraining