#!/bin/bash

#~~~~~~~~~~~~~~~ Variables to be changed by the USER ~~~~~~~~~~~~~~~~#

JOB_NAME="ALL"
REF="/scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/UCSC_sacCer3_GAP1/sacCer3.fa"
FASTQ_DIR="/scratch/cgsb/gencore/out/Gresham/2017-05-24_HKFYTBGX2/merged"
bases_trim_5prime_Read1=14
bases_trim_5prime_Read2=14
bases_trim_3prime_Read1=5
bases_trim_3prime_Read2=5
array=8-12
FASTQ_PE_1="HKFYTBGX2_n01_mini02_partii_"
FASTQ_PE_2="HKFYTBGX2_n02_mini02_partii_"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~   DO NOT EDIT BEYOND THIS LINE  ~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

RUNDIR=$PWD
DATE=`date +%m-%d-%Y`

cd $(pwd)
mkdir oe

sbatch --mail-type=ALL --mail-user=$USER@nyu.edu --job-name=$JOB_NAME -a ${array} -o $(pwd)/oe/%x_$DATE_%A_%a.out -e $(pwd)/oe/%x_$DATE_%A_%a.err cnv_detection_pipeline.sh ${REF} ${FASTQ_DIR} ${bases_trim_5prime_Read1} ${bases_trim_5prime_Read2} ${bases_trim_3prime_Read1} ${bases_trim_3prime_Read2} ${RUNDIR} ${FASTQ_PE_1} ${FASTQ_PE_2}
