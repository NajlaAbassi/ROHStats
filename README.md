# A simple step by step guide for Runs of Homozygosity analysis
## Introduction
Runs Of Homozygosity or ROH are uninterrupted homozygous regions of the genome. They appear when two copies of an ancestral haplotype are joined in an individual. Short ROH indicate distant inbreeding, while long ROH indicate recent inbreeding. This guide will help you estimate overall homozygosity in a population using PLINK data.
## What tools you will need?
+ PLINK
+ RStudio
### Step 1: Data management and Quality Control
Before you even start your analysis, you have make sure that all your datasets are in the same genome version. If not, you can use a tool like [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert your files into 1 genome assembly.
#### Quality Control
1- Remove heterosomes `plink --bfile input_data --not-chr 23 24 25 26 --make-bed --out output_data`        
In fact PLINK refers to the XX chromosome as 23, YY chromosome as 24, XY Pseudo-autosomal region of X as 25, and MT Mitochondrial as 26.          
2-Remove SNPs with MAF < 0.05, SNPs with more than 10% missing genotypes and those that divert from Hardy–Weinberg proportions with p < 0.001 `plink --bfile output_data --geno 0.1 --maf 0.05 --hwe 0.0001 --make-bed --out filtered_data`


### Contributor
+ Najla Abassi.
