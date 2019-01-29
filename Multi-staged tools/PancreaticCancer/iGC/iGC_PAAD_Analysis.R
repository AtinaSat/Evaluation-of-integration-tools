# analysis of pancreatic cancer dataset using iGC
library(iGC)
rm(list = ls())

# load data
load("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset.Rdata")
rm(paad.me)

paad.cnv[paad.cnv==2] <- 1
paad.cnv[paad.cnv==-2] <- -1
paad.cnv <- cbind(rownames(paad.cnv), paad.cnv)
colnames(paad.cnv)[1] <- "GENE"
paad.cnv[,1] <- as.character(paad.cnv[,1])

paad.ge <- cbind(rownames(paad.ge), paad.ge)
colnames(paad.ge)[1] <- "GENE"
paad.ge[,1] <- as.character(paad.ge[,1])

# Change the data types to comply with the tool requirement
library(data.table)
setDT(paad.cnv)
setDT(paad.ge)
class(paad.cnv)
class(paad.ge)

# running iGC
library(tictoc)
tic("iGC")
cna_driven_genes <- find_cna_driven_gene(paad.cnv, paad.ge)
toc()

# > cna_driven_genes <- find_cna_driven_gene(paad.cnv, paad.ge)
# Computing gain/loss sample proportion ...
# 
# Computing CNA gain driven gene records ... 
# 
# |================================| 100%
# Computing CNA loss driven gene records ... 
# 
# |================================| 100%
# > toc()
# iGC: 15.766 sec elapsed

# saving results
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/")
write.table(cna_driven_genes$gain_driven, file = "iGC_PAAD_GainDriven_genes.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$loss_driven, file = "iGC_PAAD_LossDriven_genes.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$both, file = "iGC_PAAD_BothDriven_genes.tsv", sep = "\t", quote = F,row.names = F )
