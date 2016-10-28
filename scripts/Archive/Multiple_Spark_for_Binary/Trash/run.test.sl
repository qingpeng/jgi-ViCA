#!/bin/bash

#SBATCH -p regular
#SBATCH --qos premium
#SBATCH -N 2 
#SBATCH -t 00:00:05
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
module load spark/1.5.1
start-all.sh
spark-submit --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=50g --driver-memory 50G --executor-memory 100G  /global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/Multiple_Spark_for_Binary/target/scala-2.10/simple-project_2.10-1.0.jar 
stop-all.sh

