#!/bin/bash

#SBATCH -p debug
#SBATCH -N 3
#SBATCH -C haswell
#SBATCH -t 00:30:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --image=lgerhardt/spark-2.1.1:v1
#SBATCH --volume="/global/cscratch1/sd/qpzhang>/tmpfiles:/tmp:perNodeCache=size=200G"

module load spark/2.1.1

start-all.sh

shifter spark-submit --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=120g  --driver-memory 120G --executor-memory 120G /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/spark_training_model_noPfam.py /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib  /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib.spark_model /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib.spark_scaler
shifter spark-submit --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=120g  --driver-memory 120G --executor-memory 120G /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/spark_evaluating_model.py /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_testing.svmlib  /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib.spark_model /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib.spark_scaler /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib.spark_model.report /global/projectb/scratch/qpzhang/Run_Genelearn/Full_nextflow/virus_non-virus_training.svmlib.spark_model.png


stop-all.sh
