# A simple guide for Runs of Homozygosity analysis
## Introduction
Runs Of Homozygosity or ROH are uninterrupted homozygous regions of the genome. They appear when two copies of an ancestral haplotype are joined in an individual. Short ROH indicate distant inbreeding, while long ROH indicate recent inbreeding. This guide will help you estimate overall homozygosity in a population using PLINK data.
## What tools you will need?
+ [PLINK](https://www.cog-genomics.org/plink/)
+ [RStudio](https://posit.co/download/rstudio-desktop/)

## How to use this guide?
All you need to do is clone this repository on your machine with `git clone https://github.com/NajlaAbassi/roh_analysis.git`
Then open the Rproj file and start from there.

### Step 1: Data management and Quality Control
Before you even start your analysis, you have to make sure that all your datasets are in the same genome version. If not, you can use a tool like [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert your files into the same genome assembly (hg19, hg38..).

You also need to verify the chromosomes names in the used organism. 
#### Quality Control
**1- Remove heterosomes** `plink --bfile path/to/bfile --not-chr 23 24 25 26 --make-bed --out path/to/folder`        
In fact PLINK refers to the XX chromosome as 23, YY chromosome as 24, XY Pseudo-autosomal region of X as 25, and MT Mitochondrial as 26.   
**2- Remove SNPs with MAF < 0.05, SNPs with more than 10% missing genotypes and those that divert from Hardy–Weinberg proportions with p < 0.001** `plink --bfile path/to/bfile --geno 0.1 --maf 0.05 --hwe 0.0001 --make-bed --out path/to/folder`   
### Step 2: Identification of Runs of Homozygosity
PLINK is based on the principle of a sliding window of X SNPs that traverses the entire genome to show homozygotes over a certain minimum length. THe following command can be used to call ROH:   
`--homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 2 --homozyg-gap 100 --homozyg-kb 500 --homozyg-density x --homozyg-snp y`, where x and y are respectively determined using the functions *homozyg_density()*  and *homozyg_snp()*  
### Step 3: Statistical analysis
After identifying ROH in your sample, you can now estimate homozygosity measures:   
+ Number of ROH
+ Total run length
+ FROH
+ FHOM
### Step 4: Genomic distribution of ROH (ROH islands)
ROH islands(ROHi) are regions in the genome where the percentage of individuals within a population deviates from the anticipated values predicted by a binomial distribution. To identify these regions, you can use the list of ROHs previously identified by PLINK and the [rohproc R Package](https://github.com/CeballosGene/rohproc)
### Step 5: Functional Annotation
Incorporating functional enrichment annotation allowes to explore potential connections between ROHi and specific functional pathways or gene sets. For that you can use [BiomaRt R package](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)


### Contributor
+ Najla Abassi: author and maintainer 

### Related paper
Abassi, N. et al., 2024, Investigating Inbreeding and Ancestry Dynamics in the Tunisian Population through Runs of Homozygosity Analysis ([INSERT DOI/link to preprint])
