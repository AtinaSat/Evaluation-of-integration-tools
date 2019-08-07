## Creating melanoma cancer dataset
rm(list=ls())
library(caret)

### reading copy number values and seperating them into amplification and deletion calls
cn <- read.table(gzfile("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaData/SKCM_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"), header = T, sep = "\t")
cn[1:5,1:5]
rownames(cn) <- cn[,1]
cn <- cn[,-1]

### reading gene expression file
ge.exp <- read.table(gzfile("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaData/SKCM_HiSeqV2.gz"), header = T, sep = "\t")
ge.exp[1:5,1:5]
rownames(ge.exp) <- ge.exp[,1]
ge.exp <- ge.exp[,-1]

### reading methylation data
me.exp <- read.table("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaData/SKCM.meth.by_mean.data.txt", header = T, sep = "\t")
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
skcm.ge <- t(tmp[, -rmcols])
skcm.ge <- data.frame(skcm.ge)

skcm.cnv <- cn2[rownames(skcm.ge),]
skcm.me <- me.exp2[rownames(skcm.ge),]

skcm.cnv[1:5,1:5]
skcm.ge[1:5,1:5]
skcm.me[1:5,1:5]

# changing me columns as numeric
skcm.me <- apply(skcm.me, 2, as.numeric)
rownames(skcm.me) <- row.names(skcm.ge)

# final check
all(rownames(skcm.ge) == rownames(skcm.cnv))
all(rownames(skcm.cnv) == rownames(skcm.me))

all(colnames(skcm.ge) == colnames(skcm.cnv))
all(colnames(skcm.cnv) == colnames(skcm.me))

setwd("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/")
save(file = "MelanomaRawDataset.Rdata", skcm.ge, skcm.cnv, skcm.me)
