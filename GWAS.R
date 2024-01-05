################################################################################################
######################### GWAS- GENOME WIDE ASSOCIATION STUDIES #############################

#GWAS assess the association between the trait of interest and up to millions of SNPs.
setwd("C:/Users/MARINA/Documents/MASTER/GENOMICS/ASSOCIATION STUDIES")

#---------Data and libraries ---------
library(snpStats)
ob.plink <- read.plink("obesity",)

names(ob.plink)
ob.geno <- ob.plink$genotypes
ob.geno
