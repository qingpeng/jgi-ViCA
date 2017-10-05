#!/bin/bash

#SBATCH -p debug
#SBATCH -N 2
#SBATCH -C haswell
#SBATCH -t 00:30:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --image=lgerhardt/spark-2.1.1:v1
#SBATCH --volume="/global/cscratch1/sd/qpzhang>/tmpfiles:/tmp:perNodeCache=size=200G"

module load spark/2.1.1

start-all.sh

shifter spark-submit --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=120g  --driver-memory 120G --executor-memory 120G /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/spark_training_model_noPfam.py /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/all_segment.fasta.vect.euk_and_non_euk.svmlib.0_1  /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/all_segment.fasta.vect.euk_and_non_euk.svmlib.0_1.spark_model_shift /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/all_segment.fasta.vect.euk_and_non_euk.svmlib.0_1.spark_scaler_shift

stop-all.sh
