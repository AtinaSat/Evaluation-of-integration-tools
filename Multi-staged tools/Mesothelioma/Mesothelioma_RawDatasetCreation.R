## Creating mesothelioma dataset
rm(list=ls())
library(caret)

### reading copy number values and seperating them into amplification and deletion calls
cn <- read.table(gzfile("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaData/MESO_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"), header = T, sep = "\t")
cn[1:5,1:5]
rownames(cn) <- cn[,1]
cn <- cn[,-1]

### reading gene expression file
ge.exp <- read.table(gzfile("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaData/MESO_HiSeqV2.gz"), header = T, sep = "\t")
ge.exp[1:5,1:5]
rownames(ge.exp) <- ge.exp[,1]
ge.exp <- ge.exp[,-1]

### reading methylation data
me.exp <- read.table("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaData/MESO.meth.by_mean.data.txt", header = T, sep = "\t")
me.exp[1:5,1:5]
me.exp <- me.exp[-1,]
rownames(me.exp) <- me.exp[,1]
me.exp <- me.exp[,-1]
me.exp <- na.omit(me.exp)

## identifying common genes and samples in all the datasets
a <- intersect(colnames(cn), colnames(ge.exp))
a1 <- intersect(colnames(me.exp), a)
b <- intersect(rownames(cn), rownames(ge.exp))
b1 <- intersect(rownames(me.exp), b)

ge.exp2 <- ge.exp[b1,a1]
cn2 <- cn[b1,a1]
me.exp2 <- me.exp[b1,a1]

cn2[1:5,1:5]
ge.exp2[1:5,1:5]
me.exp2[1:5,1:5]

## removing genes with zero and near zero variance in gene expression dataset and creating final dataset
# ge <-  Gene expression dataframe
# cnv <- copy number variation dataframe
# me <- methylation intensity dataframe
tmp <- t(ge.exp2)
rmcols <- nearZeroVar(tmp)
meso.ge <- t(tmp[, -rmcols])
meso.ge <- data.frame(meso.ge)

meso.cnv <- cn2[rownames(meso.ge),]
meso.me <- me.exp2[rownames(meso.ge),]

meso.cnv[1:5,1:5]
meso.ge[1:5,1:5]
meso.me[1:5,1:5]

# changing me columns as numeric
meso.me <- apply(meso.me, 2, as.numeric)
rownames(meso.me) <- row.names(meso.ge)

# final check
all(rownames(meso.ge) == rownames(meso.cnv))
all(rownames(meso.cnv) == rownames(meso.me))

all(colnames(meso.ge) == colnames(meso.cnv))
all(colnames(meso.cnv) == colnames(meso.me))

setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/")
save(file = "MesotheliomaRawDataset.Rdata", meso.ge, meso.cnv, meso.me)
