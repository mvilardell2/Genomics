BiocManager::install("dada2")
BiocManager::install("phyloseq")
BiocManager::install("ShortRead")


library(dada2)
library(ShortRead)
library(phyloseq)
library(gridExtra)
library(digest)

#Set working directory
code_path<-setwd("C:/Users/MARINA/Documents/MASTER/APPLICATIONS/Metagenomics/Session 1")
code_path
setwd(code_path)

#### Working directory to RawData 
RawDataPath <- paste(code_path,"/RawData_2k/",sep="")

#### Working directory to Databases
DDBBPath<-paste0(code_path,"/DDBB_2022_2023/")

# List fastq files  
list.files(RawDataPath)

#Sort files to ensure reads in the same order
fnFs <- sort(list.files(RawDataPath, pattern="_R1_001.fastq.gz_2k.fq.gz"))
fnRs <- sort(list.files(RawDataPath, pattern="_R2_001.fastq.gz_2k.fq.gz"))


#Infer sampe names from filenames (in which folder they are located)
# Infer sample names from filenames
sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(RawDataPath, fnFs)
fnRs <- file.path(RawDataPath, fnRs)

# Visualize the reads quality
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])


# Create a new file path to store filtered and trimmed reads
filt_path <- file.path(code_path, "DADA2/filtered") # place filtered files in filtered/subdirectory

# Define the name of output files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# FilterAndTrim-->This truncates reads at the position we specify based in the quality step
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,verbose=T)

#To see how many sequences are passing the quality filter
head(out)

#plotQualityProfile(filtFs[1:2])
#plotQualityProfile(filtRs[1:2])


#-----Errror rate estimation-----
#Each sequencing batch will have a different error rate
#The algorithm starts with the assumption that error rates are the maximum possible
# Iteratively learn an error profile from filtered reads
errF <- learnErrors(filtFs, multithread=TRUE,randomize=T,verbose=T)
errR <- learnErrors(filtRs, multithread=TRUE,randomize=T,verbose=T)

# Visualize estimated error rates 
#plotErrors(errF, nominalQ=TRUE)
#plotErrors(errR, nominalQ=TRUE)

#----------- Dereplicate fastq files to speed up computation
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample Inference
#Using the error model, we can deduce which reads are representative and which are resulting of sequencing errors
#For that, the algorithm calculates the abundance p-values for each sequence, and test the H0 that a sequence 
#with a given error rate is too abundant to be explained by sequencing errors
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Inspect the dada-class object
dadaFs[[1]]


#---------- Merge denoised reads --> it merges paired reads if they overlap exactly.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])


# ------------Organize ASVs into a sequence table (analogous to an OTU table)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


#----- Inspect distribution of sequence length of total ASV
table(nchar(getSequences(seqtab)))
#most of our reads have 439 and 440 base pairs, and really all should be like that. 
#We see that there are some sequences that are smaller or larger that the expected reads length.


# -----De novo chimera sequence detection and removal (like remove contamination)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#removeBimeraDenovo, identifies chimera sequences 

# Calculate the proportion of non-chimeric ASVs
sum(seqtab.nochim)/sum(seqtab)
#non-chimeric ASv are the filtering (ens quedem amb els non-chimeric)
#57% of non-chimeric ASV, quite low, and this is related to the fact that we are working with biopses and have a very low microbiome mass


# Calculate number of reads obtained through each step
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged",  "nonchim")
rownames(track) <- sample.names
head(track)
#for each sample we see that we started with 2000 reads, we see the filetered,....
#the non-chimeric reads will be used for the taxonomy assignment.


##---------- Assign taxonomy ##

# Path to databases
list.files(paste(code_path,"/DDBB_2022_2023/",sep=""))

# Chose taxonomy database from Analysis Configuration (Silva default)
#Taxonomy assignment to high-quality non-chimeric sequences using SILVA DATABASE
taxa.silva<-readRDS("taxa.silva.rds")
taxa.silva<- assignTaxonomy(seqtab.nochim, paste(code_path,"/DDBB_2022_2023/silva_nr_v132_train_set.fa.gz",sep=""), multithread=TRUE)
unname(head(taxa.silva))
#takes all the non-chimeric reads (high-quality) and compares to the database
#input--string of nucleotides

# Assign species if possible
taxa.silva <- addSpecies(taxa.silva, paste0(code_path,"/DDBB_2022_2023/silva_species_assignment_v132.fa.gz",""),verbose=T)
unname(head(taxa.silva))

# Replace fasta sequences from ASV name 
taxa.silva.bkp<-taxa.silva
seqtab.nochim.bkp<-seqtab.nochim

# Try also with rdp database 
#taxa.rdp<- assignTaxonomy(seqtab.nochim, paste(code_path,"/DDBB/rdp_train_set_16.fa.gz",sep=""), multithread=TRUE)
#taxa.rdp <- addSpecies(taxa.silva, paste0(code_path,"/DDBB/rdp_species_assignment_16.fa.gz",""),verbose=T)

# Phylogenetic Tree#

# Install and load packages
BiocManager::install('DECIPHER')
biocLite('DECIPHER')


library(phangorn)
library(DECIPHER)

# Run sequence alignment using DECIPHER 
sequences <- getSequences(seqtab)
names(sequences) <- sequences 
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)

## Build a neighbour-joining tree then fit a maximum likelihood tree using the neighbour-joining tree as a starting point

phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm)  # note, tip order != sequence order
fit = pml(treeNJ, data=phang_align)
# negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,
                    rearrangement='stochastic',
                    control=pml.control(trace=0))
detach('package:phangorn', unload=TRUE)

### Now we read the phylogenetic tree into a phylogenetic object (ape package)
treeSilva<-read_tree("DADA2/DADA_seq.nwk")


#################################

# Remove huge objects.
rm(derepFs)
rm(derepRs)

# Load metadata
metadata<-read.csv(paste0(code_path,"/metadata.csv"))
rownames(metadata)<-metadata$RunID

#-------------- Create PhyloSeq Object
ps_silva<-phyloseq(otu_table(seqtab.nochim.bkp, taxa_are_rows=F,
                             sample_data(metadata)),
                   tax_table(taxa.silva.bkp))
summary(ps_silva)
ps_silva

# Save the session
save.image(file=paste0(code_path,"/DADA2/DADA2_Rsession.RData"))




