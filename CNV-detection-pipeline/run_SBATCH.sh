#!/bin/bash

cd $(pwd)
mkdir oe
DATE=`date +%m-%d-%Y`
RUNDIR=$PWD
JOB_NAME="ALL"

sbatch --mail-type=ALL --mail-user=$USER@nyu.edu --job-name=$JOB_NAME -o $(pwd)/oe/%x_$DATE_%A_%a.out -e $(pwd)/oe/%x_$DATE_%A_%a.err cnv.sh /scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/UCSC_sacCer3_GAP1/sacCer3.fa /scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged 14 14 2 2 ${RUNDIR}
