### Data generation for Mesothelioma dataset using OncoDrive-CIS genes
library(reshape2)
rm(list =ls())
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset.Rdata")
rm(meso.me)

## Gene expression
topIQR <- sort(apply(meso.ge,1,IQR), decreasing = T)
topIQR2 <- topIQR[topIQR>0.9] #
meso.ge <- meso.ge[attr(topIQR2, "names"),]
meso.ge[1:5,1:5]

# Copy number data
meso.cnv <- meso.cnv[attr(topIQR2, "names"),]
meso.cnv[1:5,1:5]

all(rownames(meso.ge) == rownames(meso.cnv))
all(colnames(meso.ge) == colnames(meso.cnv))

save(meso.cnv, meso.ge, file = "/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset_OncodriveCISGenes.Rdata")
