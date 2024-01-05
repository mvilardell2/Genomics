##########################     GAS-SINGLE ASSOCIATION STUDIES      ########################################
##################################################################################################
#Genetic association studies aim to identify genetic variants that explain subject differences in qualitative or quantitative traits
#Genetic variants may include SNPs, CNVs, mosaicisms or polymorphic inversions, which can be obtained from SNP array.

setwd("~/MASTER/GENOMICS/ASSOCIATION STUDIES")
#-------Libraries----------
library(SNPassoc)

#----------- Descriptive analysis-----------
asthma<-read.delim('asthma.txt')
str(asthma,list.len=9)

#Demografic data: country, age, gender, if smoke..
#Genetic data: SNPS and Case/control variable which is the outcome variable of the analysis. 

asthma[1:5, 1:8]
#We observe that we have case-control status (0: control, 1: asthma) and another 4 variables
#encoding the country of origin, gender, age, body mass index (bmi) and smoking status (0: no smoker, 1: ex-smoker, 2: current smoker). 
#There are 51 SNPs whose genotypes are given by the alleles names.

# indicate which columns of the dataset asthma contain the SNP data
idx<-grep('^rs',colnames(asthma))
asthma.s <- setupSNP(data=asthma, colSNPs=idx, sep="")#the sep argument is to define the separator of the 2 alleles of the SNPs data.


#Descriptive analysis of SNPs
head(asthma.s$rs1422993)#information about the genotypes and the alleles.
summary(asthma.s$rs1422993)#frequency and percentatge of the genotypes and alleles and the p-value (Hardly weinberg equilibrium)
plot(asthma.s$rs1422993)

#The summary function can also be applied to the whole dataset showing the SNP labels
#with minor/major allele format, the major allele frequency the HWE test and the percentage of missing genotypes.
summary(asthma.s, print=FALSE)

#Plot missing values:
plotMissing(asthma.s, print.labels.SNPs = FALSE)
#each black line corresponds a missing in the snp.
#Missing pattern is considered random, except for the large black squares, These individuals should be cheked later. 


#--------QUALITY CONTROL----------
#WE NEED TO SEPARATE THE hwe P-VALUE IN EACH OF THE 2 GROUPS (CASE/CONTROLS)
hwe <- tableHWE(asthma.s)#HWE test--it returns the p-value
head(hwe)
#We observe that the first SNPs in the dataset are under HWE since their
#P-values rejecting the HWE hypothesis (null hypothesis) are larger than 0.05. 

#STRATIFY BY CASES ANS CONTROLS:
hwe2 <- tableHWE(asthma.s, casecontrol)#compute the p-value in each group, cases and controls

snpNHWE <- hwe2[,1]>0.05 & hwe2[,2]<0.05#Inspect the HWE in controls.Expect to find the HWE in a population that is not affected from other factors, so that is not our taget population.
rownames(hwe2)[snpNHWE]#which SNPs when separating by cases and controls reject the null hypotesis.
hwe2[snpNHWE,]
#H0: the SNP is following the HWE. 
#H1: the SNP is not following the HWE.
#rs1345267 is not in HWE within controls because its P-value is <0.05. 
#We are interested in keeping those SNPs that do not reject the null hypothesis.
#As several SNPs are tested, multiple comparisons must be considered. Threshold of 0.001

#SNPs that do not pass the HWE test must be removed form further analyses. indicate the columns of the SNPs to be kept
snps.ok <- rownames(hwe2)[hwe2[,2]>=0.001]
pos <- which(colnames(asthma)%in%snps.ok, useNames = FALSE)
asthma.s <- setupSNP(asthma, pos, sep="")


#----------------ASSOCIATION ANALYSIS ---------------------------------
#We want to find SNPs associated with asthma status encoded in the variable casecontrol.
#First association between case-control status and the fisrt SNP.
association(casecontrol ~ rs1422993, data=asthma.s )
# we observe that all genetic models but the recessive one are statistically significant.
#AIC model-->the less the better. This time, DOMINANT MODEL. 

#Conclution: this SNP is associated with asthma and, the risk of being asthmatic is 39% (OR)
# higher in perople having at least one alternative allele (T) with respect to individuals having none.
maxstat(asthma.s$casecontrol, asthma.s$rs1422993)#takes into account the 3 most common models and compute model. 

#only select one model
association(casecontrol ~ rs1422993, asthma.s, model="dominant")

#we can add some covariates for correcting the effect of the genetic variant to the case control outcome. 
association(casecontrol ~ rs1422993 + country + smoke, asthma.s)

#------------
#For multiple SNP data, our objective is to identify the variants that are significantly associated with the trait
#Fit a test for each SNP and determine which of those associations are significant. 
ans <- WGassociation(casecontrol, data=asthma.s)
head(ans)#p-value for model and for each of the snps. Keep only significant SNPS 
#perform association on those significant SNps to check results. 

#Manhattan plot of the log10 pvalue for all SNPs over all models. 
plot(ans)
#we can see the bonferroni level. 
#the overall hypotesis is whether there is any SNP that is significantly associated with the phenotype. 
#The Bonferroni correction lowers the threshold by the number of SNPs tested (0.0001=0.05/51)
#In this case, there is no SNP significant at the bonferroni level,
#so there is no SNP that is significanlty associated with asthma. 



###########################################################################################
######## EXERCISES:
#Researchers are interested in assessing possible association between candidate SNPs and the response to treatment in patients diagnosed with major depression

data <- read.delim("DM.txt")
head(data)
#Is it necessary to check HWE hypothesis for this example?
#there are not cases and controls, there are only cases, so doesn't have to to this

#Is there any SNP associated with the response to the treatment? 
#If so, write a sentence interpreting the results of that association (only for one SNP)

library(SNPassoc)
library(ggplot2) 

data.s <- setupSNP(data, 6:14, sep="")#indicate which columns correspond to the SNPs
ii <- grep("^rs|lpr", colnames(data))
data.s <- setupSNP(data, ii, sep="")

summary(data.s)

ans <- WGassociation(RESP, data.s)
ans
# get SNP name that is significant 
sig.level <- 0.05/2.2 # correct for the use of multiple genetic models
sel <- apply(ans, 1, function(x)  any(x < sig.level & !is.na(x)))
sel
names(sel[sel])

# ORs for genetic models
association(RESP ~ rs908867, data.s)
#Having at least one A allele implies 2.38 more probabilities to have a positive response to the treatment. 


##Does the result change after adjusting for other clinical covariates?
association(RESP ~ rs908867 + HDRS, data.s)
association(RESP ~ rs908867 + HDRS + PSICOT + MELANCOL + EPD_PREV, data.s)


##Create a plot with the p-values only for dominant, recessive and additive models.
data.s$RESP<-as.numeric(as.factor(data.s$RESP))
ans2 <- WGassociation(RESP, data.s, model=c("do", "re", "lo"))
plot(ans2)


  