#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=10:00:00
#SBATCH --mem=50GB

echo $(date)
################## Defining variables ###############################

#######  input variables from command line   #######
ref=$1
DATA_DIR=$2
bases_trim_5prime_Read1=$3
bases_trim_5prime_Read2=$4
bases_trim_3prime_Read1=$5
bases_trim_3prime_Read2=$6
RUNDIR=$7
file1=$8
file2=$9

####### other hard-coded variables  #######
if [ ${SLURM_ARRAY_TASK_ID} -lt 10 ]; then taskID=0${SLURM_ARRAY_TASK_ID}; else taskID=${SLURM_ARRAY_TASK_ID}; fi
ID=sample_${taskID}
PICARD_JAR='java -jar /share/apps/picard/2.8.2/picard-2.8.2.jar'
RUNDIR="$PWD"
fastq1=${DATA_DIR}/${file1}${taskID}.fastq.gz
fastq2=${DATA_DIR}/${file2}${taskID}.fastq.gz
ID=sample_${taskID}
minimum_read_length=50

############ moving in current working directory ###########
cd $RUNDIR
mkdir ${ID}
cd ${ID}

################ Data pre-processing ###################

module purge
module load trim_galore/0.4.4
module load cutadapt/intel/1.12

trim_galore --fastqc --length ${minimum_read_length} \
--trim-n \
--clip_R1 ${bases_trim_5prime_Read1} \
--clip_R2 ${bases_trim_5prime_Read2} \
--three_prime_clip_R1 ${bases_trim_3prime_Read1} \
--three_prime_clip_R2 ${bases_trim_3prime_Read2} \
-o $PWD \
--paired \
${fastq1} \
${fastq2}

module purge
module load fastqc/0.11.5
fastq1=*val_1.fq.gz
fastq2=*val_2.fq.gz
fastqc $fastq1
fastqc $fastq2

############### Alignment ####################

# performing bwa mem algorithm, sorting,indexing using samtools
module purge
module load bwa/intel/0.7.15
module load samtools/intel/1.3.1
bwa mem -t 20 $ref $fastq1 $fastq2 > ${ID}.bam
samtools sort ${ID}.bam > ${ID}.sorted.bam
samtools index ${ID}.sorted.bam

########### Alignment & insert size metrics #############
bam=${ID}.sorted.bam
# obtaining alignment metrics using Picards tools
module purge
module load picard/2.8.2
java -jar $PICARD_JAR \
CollectAlignmentSummaryMetrics \
R=$ref \
I=${bam} \
O=alignment_metrics.txt

grep -v '^#' alignment_metrics.txt | cut -f1-12 | sed '/^\s*$/d' > ${ID}_alignment_metrics.txt

# obtaining insert size metrics using Picards tools
java -jar $PICARD_JAR \
CollectInsertSizeMetrics \
INPUT=${bam}  \
OUTPUT=insert_metrics.txt \
HISTOGRAM_FILE=${ID}_insert_size_histogram.pdf

head -n 9 insert_metrics.txt | grep -v '^#' | cut -f1-19|sed '/^\s*$/d' >${ID}_insert_size_metrics.txt

############### read depth file ###################

# obtaining read depth ie coverage using samtools
module load samtools/intel/1.3.1
samtools depth -a ${bam} > ${ID}_RD.txt

############ run Pindel ###############

# obtain config file for pindel
mean_IS=$(sed -n '2p' < ${ID}_insert_size_metrics.txt | cut -f 5)
echo "${bam} ${mean_IS} ${ID}" > config_${ID}.txt

module purge
module load pindel/intel/20170402
/share/apps/pindel/20170402/intel/bin/pindel \
-T 20 \
-f $ref \
-i config_${ID}.txt \
-c ALL \
-o ${ID}_output

pindel2vcf -p ${ID}_output_D -r $ref -R UCSC_SacCer -d Feb2017 -v ${ID}_DEL_pindel.vcf
pindel2vcf -p ${ID}_output_TD -r $ref -R UCSC_SacCer -d Feb2017 -v ${ID}_DUP_pindel.vcf

module purge
module load vcftools/intel/0.1.14
vcf-concat ${ID}_DEL_pindel.vcf ${ID}_DUP_pindel.vcf > ${ID}_pindel.vcf

echo "pindel done"

############ run Lumpy ###############

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

############ run SvABA ###############

module purge
module load svaba/intel/0.2.1
svaba run -p 20 -G $ref -t $bam
cp no_id.svaba.sv.vcf ${ID}_svaba.vcf

######## Annotate all 3 VCF files #########

module load annovar/20160201
annotate_variation.pl -downdb ensGene -buildver sacCer3 sacCer3/
annotate_variation.pl --buildver sacCer3 --downdb seq sacCer3/sacCer3_seq
retrieve_seq_from_fasta.pl sacCer3/sacCer3_ensGene.txt -seqdir sacCer3/sacCer3_seq -format ensGene -outfile sacCer3/sacCer3_ensGeneMrna.fa

table_annovar.pl ${ID}_svaba.vcf sacCer3/ --buildver sacCer3 -protocol ensGene -operation g -vcfinput
table_annovar.pl ${ID}_lumpy.vcf sacCer3/ --buildver sacCer3 -protocol ensGene -operation g -vcfinput
table_annovar.pl ${ID}_pindel.vcf sacCer3/ --buildver sacCer3 -protocol ensGene -operation g -vcfinput

less ${ID}_svaba*multianno.vcf | grep -v "#" > ${ID}_SvABA_annotated.vcf
less ${ID}_lumpy*multianno.vcf | grep -v "#" > ${ID}_lumpy_annotated.vcf
less ${ID}_pindel*multianno.vcf | grep -v "#" > ${ID}_pindel_annotated.vcf

############ CNV summary and genome-wide analysis report ###############
module purge
module load r/intel/3.3.2
module load pandoc/gnu/1.19.2
cp ${RUNDIR}/CNVreport.Rmd .
Rscript -e "rmarkdown::render('CNVreport.Rmd')" ${ID}_lumpy_annotated.vcf ${ID}_pindel_annotated.vcf ${ID}_SvABA_annotated.vcf ${ID}_insert_size_metrics.txt ${ID}_alignment_metrics.txt ${ID}_RD.txt ${ID}
mv CNVreport.html ${ID}_CNVreport.html

################################################################
################## REMOVE FILES TO SAVE SPACE ##################
################################################################

# removing files to save space
ls *bam| grep -v ${ID}.sorted.bam | xargs rm
ls *bai| grep -v ${ID}.sorted.bam.bai | xargs rm

############################################################

echo "ALL DONE"
echo $(date)
