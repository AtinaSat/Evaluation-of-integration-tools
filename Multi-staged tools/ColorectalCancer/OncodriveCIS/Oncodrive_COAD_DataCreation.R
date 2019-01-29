### Data generation for colon cancer dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset.Rdata")
rm(coad.me)

## Gene expression
# Calculating interquartile range (IQR) and choosing those genes with IQR>0.55
topIQR <- sort(apply(coad.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.55]
ge.table <- coad.ge[attr(topIQR2, "names"),]
ge.table <- cbind(rownames(ge.table), ge.table)
colnames(ge.table)[1] <- "GeneID"
ge.table[1:5,1:5]


# Copy number data
coad.cnv <- coad.cnv[attr(topIQR2, "names"),]
# Changing -1 to -2 and 1 to 2 as required by the tool to include heterozygous deletion
table(coad.cnv[,2])
coad.cnv[coad.cnv==-1] <- -2
coad.cnv[coad.cnv==1] <- 2
table(coad.cnv[,2])
#changing the cnv table from wide to long
coad.cnv <- cbind(rownames(coad.cnv),coad.cnv)
colnames(coad.cnv)[1] <- "GeneID"
cnv.table <- melt(coad.cnv) 
head(cnv.table)


## Creating sample file: this file contains a list samples and a column containing 0 or 1, for normal and tumour sample, respectively
# This dataset contains no normal samples
sample.table <- data.frame(colnames(coad.cnv))
sample.table$Code <- 1
sample.table <- sample.table[-1,]
head(sample.table)

setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/OncodriveCIS/")
write.table(ge.table, file= "GeneExp_CompCOAD.tsv", sep = "\t", row.names = F, quote = F)
write.table(cnv.table, file = "CNV_CompCOAD.tsv", sep = "\t", row.names = F, quote = F)
write.table(sample.table, file = "SampleInfo_CompCOAD.tsv", sep = "\t", row.names = F, quote = F)


# Please remove column header in CNV and sample file before running Oncodrive-CIS 
# sed -i -e "1d" SampleInfo_Compcoad.tsv