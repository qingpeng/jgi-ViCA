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

shifter spark-submit --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=120g --driver-memory 120G --executor-memory 120G /global/homes/q/qpzhang/Github/jgi-ViCA/scripts/spark_prediction.py /global/projectb/scratch/qpzhang/TARA/Libsvm/all_2k.libsvm /global/homes/q/qpzhang/Github/jgi-ViCA/extra/model_tuning/Micromonas/all_segment.fasta.vect.family.all.svmlib.no4.0_1.oversampling.spark_model_shift/ /global/homes/q/qpzhang/Github/jgi-ViCA/extra/model_tuning/Micromonas/all_segment.fasta.vect.non_euk_non_virus_and_virus.svmlib.0_1.spark_scaler_shift/ /global/projectb/scratch/qpzhang/TARA/Libsvm/all_2k.libsvm.prediction_full_training

stop-all.sh
