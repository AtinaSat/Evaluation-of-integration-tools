### Data generation for pancreatic cancer dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset.Rdata")
rm(paad.me)

## Gene expression
topIQR <- sort(apply(paad.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.5]
paad.ge <- paad.ge[attr(topIQR2, "names"),]
paad.ge[1:5,1:5]

# Copy number data
paad.cnv <- paad.cnv[attr(topIQR2, "names"),]
paad.cnv[1:5,1:5]

all(rownames(paad.ge) == rownames(paad.cnv))
all(colnames(paad.ge) == colnames(paad.cnv))

save(paad.cnv, paad.ge, file = "/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset_OncodriveCISInput.Rdata")
