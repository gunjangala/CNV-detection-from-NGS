#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --mem=30GB
#SBATCH --job-name=all
#SBATCH --mail-user=ggg256@nyu.edu
#SBATCH -o /scratch/cgsb/gresham/Gunjan/oe/slurm-%A_%a.out
#SBATCH -e /scratch/cgsb/gresham/Gunjan/oe/slurm-%A_%a.err
#SBATCH --array=10-40

#############################################################################
PICARD_JAR='java -jar /share/apps/picard/2.8.2/picard-2.8.2.jar'
ref="/scratch/work/cgsb/reference_genomes/In_house/Fungi/Saccharomyces_cerevisiae/UCSC_sacCer3_GAP1/sacCer3.fa"
RUNDIR="/scratch/cgsb/gresham/Gunjan/GC_sample_analysis"
sample="sample"
ID=${sample}_${SLURM_ARRAY_TASK_ID}

cd $RUNDIR
mkdir ${ID}
cd ${ID}

cp /scratch/cgsb/gencore/out/Gresham/2016-10-26_HG7CVAFXX/merged/HG7CVAFXX_n01_mini02_${SLURM_ARRAY_TASK_ID}.fastq.gz .
cp /scratch/cgsb/gencore/out/Gresham/2016-10-26_HG7CVAFXX/merged/HG7CVAFXX_n02_mini02_${SLURM_ARRAY_TASK_ID}.fastq.gz .
gunzip *n01*
gunzip *n02*

fastq1=HG7CVAFXX_n01_mini02_${SLURM_ARRAY_TASK_ID}.fastq
fastq2=HG7CVAFXX_n02_mini02_${SLURM_ARRAY_TASK_ID}.fastq

# fastqc analysis
module purge
module load fastqc
fastqc $fastq1
fastqc $fastq2


######################### Alignment & Metrics  ###########################

# performing bwa mem algorithm, sorting,indexing using samtools
module purge
module load bwa/intel/0.7.15 
module load samtools/intel/1.3.1
bwa mem -t 16 $ref $fastq1 $fastq2 > tmp_${ID}.bam
samtools sort tmp_${ID}.bam > tmp_${ID}.sorted.bam
samtools index tmp_${ID}.sorted.bam

module purge
module load deeptools/intel/2.4.2

computeGCBias -b ${RUNDIR}/${ID}/tmp_${ID}.sorted.bam --effectiveGenomeSize 12100000 \
-g /scratch/cgsb/gresham/Gunjan/sacCer3.2bit \
-l 300 -freq ${ID}_GC.txt

correctGCBias -b ${RUNDIR}/${ID}/tmp_${ID}.sorted.bam --effectiveGenomeSize 12100000 \
-g /scratch/cgsb/gresham/Gunjan/sacCer3.2bit \
--GCbiasFrequenciesFile ${ID}_GC.txt -o ${ID}.sorted.bam

# obtaining alignment metrics using Picards tools
module purge
module load picard 
java -jar $PICARD_JAR \
CollectAlignmentSummaryMetrics \
R=$ref \
I=${ID}.sorted.bam \
O=${ID}_alignment_metrics.txt 

# obtaining insert size metrics using Picards tools
java -jar $PICARD_JAR \
CollectInsertSizeMetrics \
INPUT=${ID}.sorted.bam  \
OUTPUT=${ID}_insert_metrics.txt \
HISTOGRAM_FILE=${ID}_insert_size_histogram.pdf 

# obtaining read depth ie coverage using samtools
module load samtools/intel/1.3.1
samtools depth -a ${ID}.sorted.bam > ${ID}_RD.txt

# removing duplicates from the sorted bam file and building index using picard
java -jar $PICARD_JAR \
MarkDuplicates \
INPUT=${ID}.sorted.bam  \
OUTPUT=${ID}.rm.dup.bam \
METRICS_FILE=${ID}.rmdup.metrics.txt \
REMOVE_DUPLICATES=true

java -jar $PICARD_JAR \
BuildBamIndex \
INPUT=${ID}.rm.dup.bam

# assigning to bam(variable) the sorted bam file to be used as input while running algortihms 
bam1=${ID}.rm.dup.bam
bam=${ID}.sorted.bam

###################################################

# obtain config file for pindel
module purge
module load breakdancer/intel/1.4.5
/share/apps/breakdancer/1.4.5/intel/perl/bam2cfg.pl $bam -h > ${ID}_bd.cfg
mean_insert_size="$(less ${ID}_bd.cfg | cut -f9)"
mean_IS="$echo `expr substr $mean_insert_size 6 8`" 
echo "${bam} ${mean_IS} ${ID}" > config_${ID}.txt

module purge
module load pindel/intel/20170402
/share/apps/pindel/20170402/intel/bin/pindel \
-T 16 \
-f $ref \
-i config_${ID}.txt \
-c ALL \
-o ${ID}_output

module load python3/intel/3.5.3 
# parsing files and obtaining necessary data in separate .txt file to be read in by R
python /home/ggg256/scripts/parse_pindel_D_INV_TD_SI.py -f ${ID}_output_D > ${ID}_pindel.txt
python /home/ggg256/scripts/parse_pindel_D_INV_TD_SI.py -f ${ID}_output_TD >> ${ID}_pindel.txt
# parse pindel's output to obtain breakpoints 
python /home/ggg256/scripts/parse_pindel_novel_seq.py -f ${ID}_output_D > ${ID}_pindel_novelseq.txt
python /home/ggg256/scripts/parse_pindel_novel_seq.py -f ${ID}_output_TD >> ${ID}_pindel_novelseq.txt

pindel2vcf -p ${ID}_output_D -r $ref -R UCSC_SacCer -d Feb2017 -v ${ID}_DEL_pindel.vcf
pindel2vcf -p ${ID}_output_TD -r $ref -R UCSC_SacCer -d Feb2017 -v ${ID}_DUP_pindel.vcf

echo "pindel done"

###########################################################

# obtaining read depth ie coverage to decide bin size while running cnvnator algorithm
module load r/intel/3.3.2
Rscript /home/ggg256/scripts/read_depth_for_bin_size.R \
-r ${ID}_RD.txt \
> ${ID}_readDepth.txt
bin_size="$(cat ${ID}_readDepth.txt|replace \" ""|replace [1] "" )"

# obtaining individual fasta files from reference file in the "same directory"
# this step is very important for cnvnator to work
# also check for reference file
module load python3/intel/3.5.3 
python /home/ggg256/scripts/GAP1_fasta_to_each_chr.py 

module purge
module load cnvnator/intel/0.3.3
# predicting CNV regions
cnvnator \
-root ${ID}_out.root \
-genome $ref \
-tree $bam1 \
-unique 

# generating histogram
cnvnator \
-root ${ID}_out.root \
-genome $ref \
-tree $bam1 \
-his ${bin_size} \

# stats
cnvnator \
-root ${ID}_out.root \
-genome $ref \
-tree $bam1 \
-stat ${bin_size}

# partition
cnvnator \
-root ${ID}_out.root \
-genome $ref \
-tree $bam1 \
-partition ${bin_size}

# cnv calling
cnvnator \
-root ${ID}_out.root \
-genome $ref \
-tree $bam1 \
-call ${bin_size} > ${ID}_cnvnator.txt

module load python3/intel/3.5.3 
python /home/ggg256/scripts/parse_cnvnator.py -f ${ID}_cnvnator.txt -t deletion > ${ID}_cnv.txt
python /home/ggg256/scripts/parse_cnvnator.py -f ${ID}_cnvnator.txt -t duplication >> ${ID}_cnv.txt

# deleting the individual fasta files to save space
rm c*fa

echo "cnvnator done"

###########################################################

module purge
module load samtools/intel/1.3.1
module load bamaddrg/intel/20170402
samtools view -b -F 1294 $bam > discordants.bam
samtools view -h $bam|/share/apps/lumpy/20170331/intel/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - > splitters.bam
echo "splitters and discordants extracted"
# adding read groups in bam file(s)
bamaddrg -b $bam -s ${ID} > lumpy.merged.bam
bamaddrg -b splitters.bam -s ${ID} > lumpy.splitters.bam
bamaddrg -b discordants.bam -s ${ID} > lumpy.discordants.bam
echo "bamaddrg done"
# sort and index bam files
samtools sort lumpy.merged.bam > lumpy.sorted.bam
samtools sort lumpy.splitters.bam > splitters.sorted.bam
samtools sort lumpy.discordants.bam > discordants.sorted.bam
samtools index lumpy.sorted.bam
samtools index splitters.sorted.bam
samtools index discordants.sorted.bam
echo "lumpy sorting and indexing done"

module purge
module load lumpy/intel/20170331 
lumpyexpress \
-B lumpy.sorted.bam \
-S splitters.sorted.bam \
-D discordants.sorted.bam \
-o ${ID}_lumpy.vcf

module load python3/intel/3.5.3 
python /home/ggg256/scripts/parse_lumpy.py \
-f ${ID}_lumpy.vcf \
> ${ID}_lumpy.txt
echo "lumpy done"

###########################
mkdir SvABA
cd SvABA
module purge
module load svaba/intel/0.2.1 
svaba run -p 16 -G $ref -t ../$bam1
gunzip -c no_id.alignments.txt.gz | grep contig_name > ${ID}_plot.txt
cp no_id.svaba.sv.vcf ${ID}_svaba.vcf

module load annovar/20160201
annotate_variation.pl -downdb ensGene -buildver sacCer3 sacCer3/
annotate_variation.pl --buildver sacCer3 --downdb seq sacCer3/sacCer3_seq
retrieve_seq_from_fasta.pl sacCer3/sacCer3_ensGene.txt -seqdir sacCer3/sacCer3_seq -format ensGene -outfile sacCer3/sacCer3_ensGeneMrna.fa
table_annovar.pl ${ID}_svaba.vcf sacCer3/ --buildver sacCer3 -protocol ensGene -operation g -vcfinput
cd ${RUNDIR}/${ID}

##########################
module purge
module load r/intel/3.3.2
module load rstudio
mkdir cnv_plots_${ID}
cd cnv_plots_${ID}

cp /home/ggg256/scripts/analyzing_for_CNVs.Rmd .

Rscript -e "library(knitr); knit('analyzing_for_CNVs.Rmd')" \
-p ../${ID}_pindel.txt \
-c ../${ID}_cnv.txt \
-l ../${ID}_lumpy.txt \
-r ../${ID}_RD.txt

mv cnv.xls ${ID}_cnv.xls

######################################################################################
cd ${RUNDIR}/${ID}
mkdir chr_plots_${ID}
cd chr_plots_${ID}

cp /home/ggg256/scripts/plot_RD_whole_chr_justSample_UCSC.Rmd .
Rscript -e "library(knitr); knit('plot_RD_whole_chr_justSample_UCSC.Rmd')" \
-r ../${ID}_RD.txt \
-s ${ID}

######################################################################################
cd ${RUNDIR}/${ID}
# VOTING :
cp /home/ggg256/scripts/lab_data_voting.Rmd .

Rscript -e "library(knitr); knit('lab_data_voting.Rmd')" \
-p ${ID}_pindel.txt \
-c ${ID}_cnv.txt \
-l ${ID}_lumpy.txt \
-s ${ID}

mv lab_data_voting.md ${ID}_voting.md

######################################################################################
echo "ALL DONE"

