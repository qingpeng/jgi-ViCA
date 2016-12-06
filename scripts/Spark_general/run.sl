#!/bin/bash

#SBATCH -p debug
#SBATCH --qos normal
#SBATCH -N 2 
#SBATCH -t 00:10:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --mail-user=qingpeng@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -L projectb

module load spark
start-all.sh
spark-submit --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/ /global/homes/q/qpzhang/Bitbucket/jgi-genelearn/scripts/Spark_general/target/scala-2.11/simple-project_2.11-1.0.jar 
stop-all.sh

