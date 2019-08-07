### Data generation for Melanoma dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaRawDataset.Rdata")
rm(skcm.me)

## Gene expression
topIQR <- sort(apply(skcm.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.55] #14745 genes
ge.table <- skcm.ge[attr(topIQR2, "names"),]
ge.table <- cbind(rownames(ge.table), ge.table)
colnames(ge.table)[1] <- "GeneID"
ge.table[1:5,1:5]


# Copy number data
skcm.cnv <- skcm.cnv[attr(topIQR2, "names"),]
# Changing -1 to -2 and 1 to 2 as required by the tool to include heterozygous deletion
table(skcm.cnv[,2])
skcm.cnv[skcm.cnv==-1] <- -2
skcm.cnv[skcm.cnv==1] <- 2
table(skcm.cnv[,2])
#changing the cnv table from wide to long
skcm.cnv <- cbind(rownames(skcm.cnv),skcm.cnv)
colnames(skcm.cnv)[1] <- "GeneID"
cnv.table <- melt(skcm.cnv) 
head(cnv.table)


## Creating sample file: this file contains a list samples and a column containing 0 or 1, for normal and tumour sample, respectively
# This dataset contains no normal samples
sample.table <- data.frame(colnames(skcm.cnv))
sample.table$Code <- 1
sample.table <- sample.table[-1,]

setwd("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/OncodriveCIS/")
write.table(ge.table, file= "GeneExp_SKCM.tsv", sep = "\t", row.names = F, quote = F)
write.table(cnv.table, file = "CNV_SKCM.tsv", sep = "\t", row.names = F, quote = F)
write.table(sample.table, file = "SampleInfo_SKCM.tsv", sep = "\t", row.names = F, quote = F)


# Please remove column header in CNV and sample file before running Oncodrive-CIS 
# sed -i -e "1d" SampleInfo_SKCM.tsv