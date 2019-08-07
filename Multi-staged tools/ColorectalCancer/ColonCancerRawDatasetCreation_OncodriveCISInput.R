### Data generation for colon cancer dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset.Rdata")
rm(coad.me)

## Gene expression
# Calculating interquartile range (IQR) and choosing those genes with IQR>0.55
topIQR <- sort(apply(coad.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.55]
coad.ge <- coad.ge[attr(topIQR2, "names"),]
coad.ge[1:5,1:5]

# Copy number data
coad.cnv <- coad.cnv[attr(topIQR2, "names"),]
coad.cnv[1:5,1:5]

all(rownames(coad.ge) == rownames(coad.cnv))
all(colnames(coad.ge) == colnames(coad.cnv))

save(coad.cnv, coad.ge, file = "/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset_OncodriveCISInput.Rdata")
