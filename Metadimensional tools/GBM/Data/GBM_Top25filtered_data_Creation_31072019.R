## Code for GBM cancer data creation
## Omics data has been initially downloaded
## Reading protein data
rm(list = ls())
setwd("/home/anita/Benchmarking/GBM/Data/")

## Reading protein data
# Protein data is downloaded from https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/RPPA_RBN.gz
# These data have been normalized by RBN (replicate-base normalization) method developed by MDACC
protein.tmp <- read.delim(gzfile("GBM_RPPA_RBN.gz"), header = T, sep = "\t")
rownames(protein.tmp) <- protein.tmp[,1]
protein.tmp[,1] <- NULL
gbm.protein <- data.frame(t(protein.tmp))
gbm.protein[1:5,1:5]

## Reading mRNA array data 
# mRNA data is downloaded from https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HT_HG-U133A.gz
# This dataset shows the gene-level transcription estimates, as in log2(x+1) transformed RSEM normalized count.
rna.tmp <- read.delim(gzfile("GBM_HT_HG-U133A.gz"), header = T, sep = "\t")
rna.tmp[1:5,1:5] # 12042 obs x 540 var
rownames(rna.tmp) <- rna.tmp[,1]
rna.tmp[,1] <- NULL # 12042 obs x 539 var
gbm.rna <- data.frame(t(rna.tmp))

## Reading methylation 27k data
# Methylation data is downloaded from https://tcga.xenahubs.net/download/TCGA.GBM.sampleMap/HumanMethylation27.gz
# DNA methylation beta values are continuous variables between 0 and 1, 
# representing the ratio of the intensity of the methylated bead type to the combined locus intensity
me.tmp <- read.delim(gzfile("GBM_HumanMethylation27.gz"), header = T, sep = "\t") # 
dim(me.tmp) # 27578 obs x 289 var
rownames(me.tmp) <- me.tmp[,1]
me.tmp[1:5,1:5]
me.tmp[,1] <- NULL
me.tmp <- na.omit(me.tmp)
dim(me.tmp) # 22977 obs x 288 var
gbm.me <- data.frame(t(me.tmp))
gbm.me[1:5,1:5]

## Finding common samples in all three datasets
a <- intersect(rownames(gbm.protein), rownames(gbm.rna))
a2 <- intersect(rownames(gbm.me),a)


## Phenotypes
phenotypes <- read.delim(gzfile("GBM_clinicalMatrix.gz"), header = T, sep = "\t")
colnames(phenotypes)
head(phenotypes[,c(77, 28)])
table(phenotypes$histological_type)
table(phenotypes$GeneExp_Subtype)

rownames(phenotypes) <- gsub(pattern = "-", replacement = ".", x = phenotypes$sampleID)
gbm.pheno <- phenotypes[a2,]
table(gbm.pheno$GeneExp_Subtype)

gbm.pheno2 <- gbm.pheno

# Selecting these samples in the omic datasets
gbm.protein2 <- gbm.protein[rownames(gbm.pheno2),]
gbm.rna2 <- gbm.rna[rownames(gbm.pheno2),]
gbm.me2 <- gbm.me[rownames(gbm.pheno2),]


all(rownames(gbm.me2) == rownames(gbm.rna2))
all(rownames(gbm.rna2) == rownames(gbm.protein2))


### Data filtering
library(caret)
## Removing genes with near zero or zero variance
rm.cols <- nearZeroVar(gbm.protein2,names = TRUE) 
all.cols <- colnames(gbm.protein2)
gbm.protein2 <- gbm.protein2[,setdiff(all.cols, rm.cols)] # no proteins had zero or near zero variance

gbm.protein.var <- apply(gbm.protein2, 2, var)
top.protein.Var <- sort(gbm.protein.var, decreasing = TRUE)

## Choosing 1the top 25% of probes/genes/proteins
top25.protein <- top.protein.Var[top.protein.Var > quantile(top.protein.Var, .75)]
length(top25.protein)

par(font.axis = 2)
plot(top.protein.Var, pch = 20, ylab = "Variance", xlab = "Proteins", font.lab = 2, cex.axis =1.5, cex.lab = 1.5)
abline(h = quantile(top.protein.Var, .75), col = "blue", lwd = 3)
abline(v = length(top25.protein), col = "red", lwd = 3)

gbm.dat.protein <- gbm.protein2[,names(top25.protein)]

## Removing genes with near zero or zero variance
rm.cols <- nearZeroVar(gbm.rna2, names = T)
all.cols <- colnames(gbm.rna2)
gbm.dat.rna <- gbm.rna2[,setdiff(all.cols, rm.cols)] # no cols were removed
## Choosing 1the top 25% of probes/genes/proteins
mrna.Var <- apply(gbm.dat.rna, 2, var)
top.mrna.Var <- sort(mrna.Var, decreasing = T)

top25.mrna <- top.mrna.Var[top.mrna.Var > quantile(top.mrna.Var, .75)]

par(font.axis = 2)
plot(top.mrna.Var, pch = 20, ylab = 'Variance', xlab = "mRNA genes", font.lab = 2, cex.axis =1.5, cex.lab = 1.5)
abline(h = quantile(top.mrna.Var, .75), col = "blue", lwd = 3)
abline(v = length(top25.mrna), col = "red", lwd = 3)

gbm.dat.rna <- gbm.dat.rna[,attr(top25.mrna, "names")]

## Removing probes with near zero or zero variance
rm.cols <- nearZeroVar(gbm.me2, names = T) # no cols were removed
all.cols <- colnames(gbm.me2)
gbm.dat.me <- gbm.me2[,setdiff(all.cols, rm.cols)]
## Choosing the top 25% of probes/genes/proteins
me.Var <- apply(gbm.dat.me, 2, var)
top.me <- sort(me.Var, decreasing = T)
top25.me <- top.me[top.me > quantile(top.me, .75)]

par(font.axis = 2)
plot(top.me, pch = 20, ylab = "Variance",  xlab = "Methylation probes", font.lab = 2, cex.axis =1.5, cex.lab = 1.5)
abline(h = quantile(top.me, .75), col = "blue", lwd = 3)
abline(v = length(top25.me), col = "red", lwd = 3)

gbm.dat.me <- gbm.dat.me[,attr(top25.me, "names")]

#########
dim(gbm.dat.protein)
# [1] 100  33
dim(gbm.dat.rna)
# [1]  100 3011
dim(gbm.dat.me)
# [1]  100 5744

# Creating a file that contains the samples and labels (Classical, Mesenchymal, Neural and Proneural)
# The file contains two columns: Samples and Labels
data.labels <- gbm.pheno2[,c(1,28)]
head(data.labels)
colnames(data.labels) <- c("Sample", "Label")
data.labels$Sample <- gsub(pattern = "-", replacement = ".", x = data.labels$Sample)
data.labels$Label <- as.character(data.labels$Label)

head(data.labels)
table(data.labels$Label)
# Classical Mesenchymal      Neural   Proneural 
#   28          24            17          31 

setwd("/home/anita/Benchmarking/GBM/Data/")
write.table(x = data.labels, file = "GBMSamplesAndLabels.txt", sep ="\t",
            row.names = F, quote = F)

# Saving the data for future use
setwd("/home/anita/Benchmarking/GBM/Data/")
saveRDS(gbm.dat.me, file = "GBM_Methylation_Top25Filtered.rd", version = 2)
saveRDS(gbm.dat.protein, file = "GBM_Protein_Top25Filtered.rd", version = 2)
saveRDS(gbm.dat.rna, file = "GBM_RNA_Top25Filtered.rd", version = 2)
saveRDS(data.labels, file = "GBM_DataLabels.rd", version = 2)
save(gbm.dat.rna,gbm.dat.protein,gbm.dat.me,data.labels, file = "GBM_Top25filtered_data.RData")
