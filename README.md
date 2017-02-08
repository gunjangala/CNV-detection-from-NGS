### Detecting copy-number variants present in a genome (here,yeast) from next generation sequencing data

Copy number variants (CNVs) ie duplications and deletions of genome sequences encompass genes leading to dosage imbalances resulting in phenotypic diversity and rapid adaptive evolution. Thus, precise characterization of novel sequences resulting in CNVs also known as breakpoints is very crucial for assessing their functional impact. The strong demand for next generation sequencing based CNV analysis has fuelled development of numerous computational methods and tools for CNV detection. CNVs typically result in two novel features in DNA sequence that are exploited by detection approaches: 1) an increase, or decrease, in the total amount of unique sequence and 2) generation of novel sequence that is not present in the genome and therefore not mapped using reference-based alignment. We have used algorithms based on three complementary strategies for detecting CNVs using short read DNA sequencing termed read depth (RD) mapping, discordant read (DR) pair mapping, and split read (SR) mapping to detect breakpoints in microbial populations propagated in chemostats under nutrient-limited conditions. Prior to its application on this research data, we studied the effect of different read depth on performance of these algorithms on simulated data and found that reducing read-depth reduces the sensitivity of all algorithms. As these algorithms are based on distinct different approaches, we observed that integrating all three algorithms in a computational pipeline improved sensitivity at a lower read depth of 10X. Implementing this pipeline on research data helps us explore the diversity of different CNV alleles and infer the mechanisms underlying their formation and occurrence in evolved populations.
(Algorithms used: Lumpy(SR+RP), Cnvnator(RD) and Pindel(SR))

These scripts can be used for following purpose(s) :

**Simulations:**
  + simulate and detect structural variants(SV) in reference genomes (.pbs file), 
  + evaluate performances of tools/packages used in SV detection based on False discovery rate (FDR), precision, sensitivity and F-score (.Rmd file), and
  + gives an overview of how many tools successfully detects total of simulated SVs.
  
**Real/Lab data:**
  + detect SVs present in yeast genomes.
  + plot read depth vs position for variants detected in genomes as well as whole chromosomes.
  + gives an overview of how many tools successfully detects total of SVs present in genomes.
