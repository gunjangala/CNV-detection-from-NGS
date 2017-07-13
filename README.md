# An integrated approach to detect copy number variants in next-generation sequencing data

Copy number variants (CNVs) encompass genes leading to dosage imbalances resulting in phenotypic diversity and rapid adaptive evolution. Numerous computational methods have been developed for CNV detection. CNVs typically result in two novel features in DNA sequence that are exploited by detection approaches: 1) an increase, or decrease, in the total amount of unique sequence and 2) generation of novel sequence not present in the reference genome and therefore not mapped using reference-based alignment. We tested algorithms based on three complementary strategies for detecting CNVs using short read DNA sequencing: read depth (RD) mapping, discordant read (DR) pair mapping, and split read (SR) mapping. We studied the effect of read depth on algorithm performance using simulated data and found that reducing read-depth reduces the sensitivity of all algorithms. However, as different algorithms are based on different approaches, integrating results from multiple algorithms in a computational pipeline improved sensitivity compared to any single algorithm. We used our integrated pipeline to detected CNVs that emerged in microbial populations propagated in chemostats under nutrient-limited conditions. Implementing this pipeline provides a means of exploring the diversity of different CNV alleles and inferring the mechanisms underlying their formation in evolved populations.

Algorithms used: [Lumpy(SR/RP)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84),
[Pindel(SR)](https://academic.oup.com/bioinformatics/article/25/21/2865/2112044/Pindel-a-pattern-growth-approach-to-detect-break)and 
[SvABA](http://biorxiv.org/content/early/2017/02/01/105080)


Short information regarding contents of sub-folders(s):

```other``` (simulation study):
    - simulate and detect CNVs in reference genome
    - evaluate performances of tools/packages used in SV detection based on False discovery rate (FDR), precision, sensitivity and F-score
    - some old .pbs scripts

```pipeline``` (pipeline implemented on genomics data):
    - detect CNVs present in genomes,
    - plot read depth vs position for portions of genome as well as whole chromosomes, and
    - gives an overview of how many tools successfully detect total of CNVs present in genomes.

This pipeline is tested and optimized on haploid yeast genome (so far). 
To run this pipeline, we require:
  - paired-end FASTQ files as input, and
  - all files from github folder copied to your working/project directory 

## CNV detection pipeline

### Files description:

  - ```run_SBATCH.sh``` : This is a wrapper script. Here, one should change the following parameters before executing the script:
    - job name
    - Path to reference directory
    - Path to directory containing paired-end FASTQ files
    - number of bases to trim from 5' end of Read1 and Read2
    - number of bases to trim from 3' end of Read1 and Read2
    - range of array jobs to run (eg. 1-10). This number should be in accordance with nomenclature adapted by fastq files.
    - FASTQ file name prefix for each read pair. For example: 
    
      If paired-end file names are, ```HKFYTBGX2_n01_mini02_partii_10.fastq.gz``` and ```HKFYTBGX2_n02_mini02_partii_10.fastq.gz```, set the variables 
        - FASTQ_PE_1 as **"HKFYTBGX2_n01_mini02_partii_"** and 
        - FASTQ_PE_2 as **"HKFYTBGX2_n02_mini02_partii_"**
       
  Example code (for setting variables):
  
  ```
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
 ```
    
- ```cnv_detection_pipeline.sh``` : This script contains steps of actual pipeline. Do not edit this script. 

- ```CNVreport.Rmd``` : This Rmarkdown file does all downstream analysis. It generates a HTML output for each sample processed.

- ```consolidate_all_samples_results.sh``` : Run this as follows if you want to look at CNV output chromosome-wise (so far adapted for just chr4,8 and 11) for each tool used. This step is **not** mandatory/required.
``` 
cd <name of project directory>
sh consolidate_all_samples_results.sh
```

### Usage

Create a new directory on HPC Prince cluster and enter in this directory as follows:

```
mkdir <name of project directory>
cd <name of project directory>
```

Run the following code on HPC Prince to download the pipeline. Run the following code to get started.
```
git clone https://github.com/gunjangala/CNV-detection-from-NGS
cp CNV-detection-from-NGS/pipeline/* .
rm -rf CNV-detection-from-NGS/
```

Set correct variables in ```run_SBATCH.sh``` and then run the pipeline as follows:
```
sh run_SBATCH.sh
```

Check for successful submission and running of jobs as follows:
```
watch squeue -u <your net ID>
``` 
To understand what going on, refer [this](https://wikis.nyu.edu/display/NYUHPC/Slurm+Tutorial) website and go to ```Check job status``` tab.

### HTML output
If you wish to extract just HTML output of all samples, run the following code:
```
cd <name of project directory>
mkdir HTML_OUTPUTS
cp sample*/*_CNVreport.html HTML_OUTPUTS/.
```

```HTML_OUTPUTS``` will now contain output HTML files for all samples processed.

### O/E logs:
All the ERROR and OUTPUT logs for the run could be found in ```<name of project directory>/oe``` directory
Eg.

    - error log file: ALL_1372744_10.err
    - output log file: ALL_1372744_10.out
    
Here, 

    - "ALL" refers to the job name
    - "1372744" refers to job ID
    - "10" refers to array ID and/or sample number
    - file extension ie "err" stands for error log and "out" stands for output log.
    

***
