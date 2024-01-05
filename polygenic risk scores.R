##########################     POLYGENIC RISK SCORES      ########################################
##################################################################################################
setwd("~/MASTER/GENOMICS/ASSOCIATION STUDIES")

#------------------1. GENETIC SCORES------------------------
#GENETIC risk scores is given by the number of risk alleles that an individual carries.
#They aim to assess the collective prediction of the phenotype risk associated with mutliple SNPs.

#------- Load libraries-----------
library(SNPassoc)
library(MASS)
library(PredictABEL)
library(ggplot2)

#------- Association analysis ------------
asthma <- read.delim("asthma.txt")
idx <- grep("^rs", colnames(asthma))
asthma.s <- setupSNP(data=asthma, colSNPs=idx, sep="")
ans <- WGassociation(casecontrol, asthma.s, model="log-add")

# We choose the SNPs that show an association at 0.1 level to be included in a multivariate model:
sel <- labels(asthma.s)[additive(ans)<0.1]
asthma.sel <- asthma.s[,sel]
head(asthma.sel)

#We join the SNPs and case/control variable in a single data.frame. Using the additive(), we recode the SNP
#genotypes as 0, 1 or 2 where the homozygous for the alternative allele is 2.
#As such, the genetic score is the sum of alternative alleles over selected SNPs

asthma.sel <- data.frame(lapply(asthma.sel, additive)) #lapply takes each column of the data frame, and it applies the additive function--it returns a list--so we transform to a data.frame
dd.end <- data.frame(casecontrol=asthma.s$casecontrol, asthma.sel) #join the SNPs scores and the case/control variable
head(dd.end)

#A multivariate regression model is fitted with all the SNPs to further select those that best predict the asthma status.
#An automatic variable selection, based on Akaike information criteria (AIC), can be applied to the model using the
#stepAIC() function from library MASS. Note that the function does not support missing values and requires having complete cases (i.e. no missing values on genotypes).

library(MASS)
dd.end.complete <- dd.end[complete.cases(dd.end),] #select only those individuals with all information (we don't want missings)
mod <- stepAIC(glm(casecontrol ~ ., dd.end.complete,
                   family="binomial"),
               method="forward", trace=0)


#---------- Genetic score computation ---------------

#The selected SNPs, from the multivariate model, in snps.score are then used to create the genetic score. 
#First, a summary shows the SNPs’ associations of the selected SNPs with the trait

summary(mod)#we are selecting the SNPs that best predict the asthma status and which is their effect. (so we are removing the other SNPs)

#Get the names of the selected SNPs (the important ones that explain the asthma status)
snps.score <- names(coef(mod))[-1]#with [-1] we are removing the intercept column because is of no interest
snps.score

#The function riskScore() of the PredictABEL package creates the risk score:
pos <- which(names(dd.end.complete)%in%snps.score)
pos#returns the position of the names in dd.end.complete that are in snps.score

score <- riskScore(mod, data=dd.end.complete, 
                   cGenPreds=pos,
                   Type="unweighted")
table(score)
#in this case, we are not considering the effects of the alleles (mod)

barplot(table(score))

#--------- Association analysis includiing genetic Scores ----------

#Once the genetic score is created, we test its association with asthma status using a general linear model:
mod.lin <- glm(casecontrol~score, dd.end.complete,
               family="binomial")
mod.lin
exp(mod.lin$coefficients[2])#We observe that the risk of asthma increases 24% per each risk allele.

#The predictive power of the genetic score can be assessed by computing the area under the ROC curve (AUC).
predrisk <- predRisk(mod.lin, dd.end.complete)
plotROC(data=dd.end.complete, cOutcome=1,
        predrisk = predrisk)

##################################################################################
## EXERCISE
#Researchers are now interested in creating a genetic score to predict the response to 
#treatment in patients diagnosed with major depression (file ‘DM.txt’ - NOTE: check how alleles are separated).

#------- Libraries -----------------
library(SNPassoc)
library(MASS)
library(PredictABEL)
library(ggplot2)

#-----------Association analysis---------------
#Get data
DM<-read.delim('DM.txt')
head(DM)
#Indicate which are the SNPs columns
DM.s<-setupSNP(data=DM,colSNPs = 6:14, sep="")#we are saying to R which are the SNP columns

#Association analysis
ans<-WGassociation(RESP,DM.s, model='log-add')

##--Select SNPs
sel<-labels(DM.s)[additive(ans)<0.1]#to follow recomendations
sel<-labels(DM.s)#in this case this is better since we are dealing with few SNPs
DM.sel<-DM.s[,sel]#data frame with the columns of SNPs selected
head(DM.sel)

##--Genotype to 0,1 or 2
data.sel <- data.frame(lapply(DM.sel, additive))
dd.end <- data.frame(response=DM.s$RESP, data.sel)
head(dd.end)
table(dd.end$response)

##Generate a model to select those SNPs that better predict the response to treatment.

dd.end.complete <- dd.end[complete.cases(dd.end),]#Just not select the missing values patients

dd.end.complete$response<-as.numeric(as.factor(dd.end.complete$response))-1#we need the variable treatment to be numeric
table(dd.end.complete$response)

mod <- stepAIC(glm(response ~ ., dd.end.complete,
                   family="binomial"),
               method="forward", trace=0)

#-------------------- GENETIC SCORE COMPUTATION ---------------------
summary(mod)

#We select the SNPs that is important to create the genetic score: 
snps.name<-names(coef(mod))[-1]

#Position of the columns of SNP selected: 
pos<-which(names(dd.end.complete)%in%snps.name)
pos

#Risk score: 
score <- riskScore(mod, data=dd.end.complete, 
                   cGenPreds=pos,
                   Type="unweighted")
table(score)


#-----genetic scores
mod.lin <- glm(response~score, dd.end.complete,
               family="binomial")
summary(mod.lin)
exp(coef(mod.lin)[2])#the probability of a positive response to treatment increases 94% per each allele.  

#For instance, the risk of people having 4 risk alleles vs 2 risk alleles:
exp(2*coef(mod.lin)[2])


predrisk <- predRisk(mod.lin, dd.end.complete)
plotROC(data=dd.end.complete, cOutcome=1,
        predrisk = predrisk)
