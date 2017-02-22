### Detecting copy-number variants present in a genome (here,yeast) from next generation sequencing data

Copy number variants (CNVs) encompass genes leading to dosage imbalances resulting in phenotypic diversity and rapid adaptive evolution. Numerous computational methods have been developed for CNV detection. CNVs typically result in two novel features in DNA sequence that are exploited by detection approaches: 1) an increase, or decrease, in the total amount of unique sequence and 2) generation of novel sequence not present in the reference genome and therefore not mapped using reference-based alignment. We tested algorithms based on three complementary strategies for detecting CNVs using short read DNA sequencing: read depth (RD) mapping, discordant read (DR) pair mapping, and split read (SR) mapping. We studied the effect of read depth on algorithm performance using simulated data and found that reducing read-depth reduces the sensitivity of all algorithms. However, as different algorithms are based on different approaches, integrating results from multiple algorithms in a computational pipeline improved sensitivity compared to any single algorithm. We used our integrated pipeline to detected CNVs that emerged in microbial populations propagated in chemostats under nutrient-limited conditions. Implementing this pipeline provides a means of exploring the diversity of different CNV alleles and inferring the mechanisms underlying their formation in evolved populations.

Scripts in this folder can be used for following purpose(s):

(Algorithms used: [Lumpy(SR+RP)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84)
Lumpy(SR+RP), Cnvnator(RD), Pindel(SR) and SvABA)

**Simulations:**
  + simulate and detect structural variants(SV) in reference genomes (.pbs file), 
  + evaluate performances of tools/packages used in SV detection based on False discovery rate (FDR), precision, sensitivity and F-score (Rscript),
  + gives an overview of how many tools successfully detects total of simulated SVs (Rscript), and
  + Analysis for cnv detection in mixed population samples (Rscript).
  
**Real/Lab data:**
  + detect CNVs present in yeast genomes (.pbs file),
  + plot read depth vs position for variants detected in genomes as well as whole chromosomes (Rscript), and
  + gives an overview of how many tools successfully detects total of SVs present in genomes (Rscript).
