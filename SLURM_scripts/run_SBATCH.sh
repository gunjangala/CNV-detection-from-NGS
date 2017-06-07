#!/bin/bash

cd $(pwd)
mkdir oe
JOB_NAME="ALL"
DATE=`date +%m-%d-%Y`

sbatch --mail-type=ALL --mail-user=$USER@nyu.edu --job-name=$JOB_NAME -o $(pwd)/oe/%x_$DATE_%A.out -e $(pwd)/oe/%x_$DATE_%A.err /home/ggg256/scripts/5JUN_01to09.sh $(pwd)

echo $DATE
echo $JOB_NAME

