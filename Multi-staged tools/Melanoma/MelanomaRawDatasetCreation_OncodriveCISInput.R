### Data generation for Melanoma dataset for OncoDrive-CIS analysis
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaRawDataset.Rdata")
rm(skcm.me)

## Gene expression
topIQR <- sort(apply(skcm.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.55] #14745 genes
skcm.ge <- skcm.ge[attr(topIQR2, "names"),]
skcm.ge[1:5,1:5]

# Copy number data
skcm.cnv <- skcm.cnv[attr(topIQR2, "names"),]
skcm.cnv[1:5,1:5]

all(rownames(skcm.ge) == rownames(skcm.cnv))
all(colnames(skcm.ge) == colnames(skcm.cnv))

save(skcm.cnv, skcm.ge, file = "/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaRawDataset_OncodriveCISInput.Rdata")

