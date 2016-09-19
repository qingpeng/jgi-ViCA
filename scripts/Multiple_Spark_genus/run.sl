#!/bin/bash

#SBATCH -p regular
#SBATCH --qos premium
#SBATCH -N 20 
#SBATCH -t 20:00:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --mail-user=qingpeng@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --ccm
module load spark
start-all.sh
spark-submit --master $SPARKURL --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/  --conf spark.driver.maxResultSize=50g --driver-memory 50G --executor-memory 100G /global/u2/q/qpzhang/Dropbox/Bitbucket/jgi-genelearn/scripts/Multiple_Spark_genus/target/scala-2.10/simple-project_2.10-1.0.jar 
stop-all.sh

