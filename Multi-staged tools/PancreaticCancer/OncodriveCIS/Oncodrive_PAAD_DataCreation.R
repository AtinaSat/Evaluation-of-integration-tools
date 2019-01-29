### Data generation for pancreatic cancer dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("PancreaticCancerRawDataset.Rdata")
rm(paad.me)

## Gene expression
topIQR <- sort(apply(paad.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.5]
ge.table <- paad.ge[attr(topIQR2, "names"),]
ge.table <- cbind(rownames(ge.table), ge.table)
colnames(ge.table)[1] <- "GeneID"
ge.table[1:5,1:5]


# Copy number data
paad.cnv <- paad.cnv[attr(topIQR2, "names"),]
# Changing -1 to -2 and 1 to 2 as required by the tool to include heterozygous deletion
table(paad.cnv[,2])
paad.cnv[paad.cnv==-1] <- -2
paad.cnv[paad.cnv==1] <- 2
table(paad.cnv[,2])
#changing the cnv table from wide to long
paad.cnv <- cbind(rownames(paad.cnv),paad.cnv)
colnames(paad.cnv)[1] <- "GeneID"
cnv.table <- melt(paad.cnv) 
head(cnv.table)


## Creating sample file: this file contains a list samples and a column containing 0 or 1, for normal and tumour sample, respectively
# This dataset contains no normal samples
sample.table <- data.frame(colnames(paad.cnv))
sample.table$Code <- 1
sample.table <- sample.table[-1,]

setwd("PancreaticCancer/OncodriveCIS/")
write.table(ge.table, file= "GeneExp_CompPAAD.tsv", sep = "\t", row.names = F, quote = F)
write.table(cnv.table, file = "CNV_CompPAAD.tsv", sep = "\t", row.names = F, quote = F)
write.table(sample.table, file = "SampleInfo_CompPAAD.tsv", sep = "\t", row.names = F, quote = F)


# Please remove column header in CNV and sample file before running Oncodrive-CIS 
# sed -i -e "1d" SampleInfo_CompPAAD.tsv
