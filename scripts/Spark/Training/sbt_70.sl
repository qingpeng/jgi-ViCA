#!/bin/bash

#SBATCH -p regular
#SBATCH --qos normal 
#SBATCH -N 20 
#SBATCH -t 35:00:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --ccm
module load spark
start-all.sh
spark-submit --master $SPARKURL --conf spark.eventLog.enabled=true --conf spark.eventLog.dir=$SCRATCH/spark/ --driver-memory 15G --executor-memory 32G /global/homes/q/qpzhang/Dropbox/NewBitbucket_for_GeneLearn/jgi-genelearn/scripts/Spark/Training/target/scala-2.10/simple-project_2.10-1.0.jar /global/projectb/scratch/qpzhang/Full_Training/Pfam/training.vect.all.sort.log /global/projectb/scratch/qpzhang/Full_Training/Pfam/testing.vect.all.sort.log /global/projectb/scratch/qpzhang/Full_Training/Pfam/Test_GBT_Scaling_70/ sbt 70  
stop-all.sh

