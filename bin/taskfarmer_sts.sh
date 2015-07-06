#!/bin/sh

. env.sh
module load biopython
#Designed to be submitted to 8 core commodity nodes with two threads per worker
for i in {1..4}
do
   echo "start worker $i"
   tfmq-worker -z -q singletaxontraining &
done
wait