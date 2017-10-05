#!/bin/bash

#SBATCH -p debug
#SBATCH -N 2
#SBATCH -C haswell
#SBATCH -t 00:30:00
#SBATCH -e mysparkjob_%j.err
#SBATCH -o mysparkjob_%j.out
#SBATCH --image=lgerhardt/spark-2.1.1:v1
#SBATCH --volume="/global/cscratch1/sd/<user_name>/tmpfiles:/tmp:perNodeCache=size=200G"

module load spark/2.1.1

start-all.sh

shifter spark-submit $SPARK_EXAMPLES/python/pi.py

stop-all.sh
