## Code for hepatocellular carcinoma data creation
# we select top 10% for miRNA, mRNA and top 5% of methylation
## Omics data has been initially downloaded
## Reading miRNA data
rm(list = ls())
setwd("/home/anita/Benchmarking/HepatocellularCarcinoma/data/GSE76903_RAW_miRNA/")
file_list <- list.files()
dat.tmp <- do.call("cbind",lapply(file_list, FUN=function(files){read.table(gzfile(files),header=FALSE, sep="\t")}))
dat.tmp <- dat.tmp[-1,]
colno <- seq(2,120,2)
dat.mirna <- dat.tmp[,colno]
rownames(dat.mirna) <- dat.tmp[,1]
col.names <- substr(file_list, 12,14)
col.names <- gsub("[.]", "", col.names)
colnames(dat.mirna) <- col.names
dat.mirna <- t(dat.mirna)
dat.mirna <- apply(dat.mirna, 2, as.numeric)
rownames(dat.mirna) <- col.names
dat.mirna <- data.frame(dat.mirna)
class(dat.mirna[,1])
dat.mirna[1:5,1:5]

## Reading mRNA data
setwd("../GSE77509_RAW_mRNA/")
file_list2 <- list.files()
dat.tmp <- do.call(cbind,lapply(file_list2, FUN=function(files){read.table(gzfile(files), header = FALSE, sep = "\t")}))
dat.tmp <- dat.tmp[-1,]
dat.mrna <- dat.tmp[,colno]
rownames(dat.mrna) <- dat.tmp[,1]
col.names <- substr(file_list2, 12, 14)
col.names <- gsub("[.]", "", col.names)
colnames(dat.mrna) <- col.names
dat.mrna <- t(dat.mrna)
dat.mrna <- apply(dat.mrna, 2, as.numeric)
rownames(dat.mrna) <- col.names
dat.mrna <- data.frame(dat.mrna)
dat.mrna[1:5,1:5]

## Reading methylation data
setwd("../GSE77269_RAW_methylation/")
file_listtmp <- list.files()
file_list3 <- file_listtmp[grep("txt.",file_listtmp)]
dat.tmp <- do.call(cbind,lapply(file_list3, FUN=function(files){read.table(gzfile(files), header = FALSE, sep = "\t")}))
colno <- seq(2,180,3)
dat.me <- dat.tmp[,colno]
rownames(dat.me) <- dat.tmp[,1]
dat.me <- na.omit(dat.me)
col.names <- substr(file_list3, 12, 14)
col.names <- gsub("[.]", "", col.names)
colnames(dat.me) <- col.names
dat.me <- data.frame(t(dat.me))
dat.me[1:5,1:5]

## Checking if the sample names are same in all the datasets
all(rownames(dat.me) == rownames(dat.mirna))
all(rownames(dat.mrna) == rownames(dat.mirna))

# log transforming the data (only miRNA and mRNA) prior to selection of features based on variance
dat.mrna <- log(dat.mrna+1)
dat.mirna <- log(dat.mirna+1)

dat.mirna = dat.mirna[,colSums(dat.mirna==0) < 60*.25]
dat.mrna = dat.mrna[,colSums(dat.mrna==0) < 60*0.25]
dat.me = dat.me[,colSums(dat.me==0) < 60*0.25]

### Data filtering
library(caret)
## Removing genes with near zero or zero variance
rm.cols <- nearZeroVar(dat.mirna,names = TRUE)
all.cols <- colnames(dat.mirna)
mirna <- dat.mirna[,setdiff(all.cols, rm.cols)]
## Selection of most variant genes and probes
mirna.Var <- apply(mirna, 2, var)
top.mirna.Var <- sort(mirna.Var, decreasing = TRUE)
#plots
par(font.axis = 2)
plot(top.mirna.Var, pch = 20, ylab = 'Variance', xlab = "miRNA genes", font.lab = 2, cex.axis =1.5, cex.lab = 1.5)
abline(h = quantile(mirna.Var, .90), col = "blue", lwd = 3)
abline(v = length(top.mirna.Var[top.mirna.Var > quantile(top.mirna.Var, .90)]), col = "red", lwd = 3)
# Selecting top 10%
top10.mirna <-  top.mirna.Var[top.mirna.Var > quantile(top.mirna.Var, .90)]
mirna <- mirna[,names(top10.mirna)]

## Removing genes with near zero or zero variance
rm.cols <- nearZeroVar(dat.mrna, names = T)
all.cols <- colnames(dat.mrna)
mrna <- dat.mrna[,setdiff(all.cols, rm.cols)]
mrna.Var <- apply(mrna, 2, var)
top.mrna.Var <- sort(mrna.Var, decreasing = T)
#plots
par(font.axis = 2)
plot(top.mrna.Var, pch = 20, ylab = 'Variance', xlab = "mRNA genes", font.lab = 2, cex.axis =1.5, cex.lab = 1.5)
abline(h = quantile(top.mrna.Var, .90), col = "blue", lwd = 3)
abline(v = length(top.mrna.Var[top.mrna.Var > quantile(top.mrna.Var, .90)]), col = "red", lwd = 3)
# Selection top 10%
top10.mrna <-  top.mrna.Var[top.mrna.Var > quantile(top.mrna.Var, .90)]
mrna <- mrna[,names(top10.mrna)]

## Top quartile of most variant genes and probes
rm.cols <- nearZeroVar(dat.me, names = T)
all.cols <- colnames(dat.me)
me <- dat.me[,setdiff(all.cols, rm.cols)]
#me <- dat.me
me.Var <- apply(me, 2, var)
top.me <- sort(me.Var, decreasing = T)
#plots
par(font.axis = 2)
plot(top.me, pch = 20, ylab = 'Variance', xlab = "Methylation probes", font.lab = 2, cex.axis =1.5, cex.lab = 1.5)
abline(h = quantile(top.me, .95), col = "blue", lwd = 3)
abline(v = length(top.me[top.me > quantile(top.me, .95)]), col = "red", lwd = 3)
# Taking top 5%
top5.me <- top.me[top.me>quantile(me.Var, .95)]
me <- me[,names(top5.me)]

#########
setwd("../")
dim(mirna)
dim(mrna)
dim(me)
# Reading the file that contains the labels (Normal, tumour or PVTT). The file contains two columns: Samples and Labels
data.labels <- read.table("SamplesAndLabels.txt", header = T, sep ="\t")

# Saving the data for future use
save(me,mirna,mrna, data.labels, file = "HPC_Top10Top5Filtered_data.Rdata")


