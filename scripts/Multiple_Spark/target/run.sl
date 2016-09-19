#!/bin/bash

#SBATCH -p regular
#SBATCH --qos premium
#SBATCH -N 40 
#SBATCH -t 04:00:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --ccm
module load spark
start-all.sh
spark-submit --master $SPARKURL --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=20g --driver-memory 25G --executor-memory 30G /global/u2/q/qpzhang/Dropbox/Bitbucket/jgi-genelearn/scripts/Multiple_Spark/target/scala-2.10/simple-project_2.10-1.0.jar 
stop-all.sh

