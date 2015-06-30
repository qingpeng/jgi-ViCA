#!/bin/bash -l

# usage genemarksub.sh <1 abs input dir path> <2 abs running dir path> <3 params> <4 file number for testing>
# Adam Rivers US DOE JGI
#$ -w e

#source genemark program
source env.sh

#set human friendly variable names
inputdir="$1"
echo "input directory set as " "$inputdir"
runningdir="$2"
echo "runningdir set as " "$runningdir"
params="$3"
echo "params set as " $params
cd $runningdir

# look for optional test id or else use the UGE TASK ID 
taskid=$4
if [ -z "$taskid" ]
then
	taskid="$SGE_TASK_ID"
fi
if [ -z "$taskid" ]
then
	echo "task id not defined" >&2
	exit 1
fi

#create and move into a temporary directory for the specific genome
mkdir temp"$taskid"
cd temp"$taskid"

#convert the taskid into the name of a file to be run
filename=`ls -1 $inputdir | tail -n +"$taskid" | head -1`

#Run command: genemark + parameters from create_training_data.py + thread number
$params "--output" "$filename"  "--name" "$filename" "$inputdir/$filename" 
# move the model file to the correct directory
mv *hmm.mod $runningdir/.

#clean up directory
cd ..
rm -r temp"$taskid"

