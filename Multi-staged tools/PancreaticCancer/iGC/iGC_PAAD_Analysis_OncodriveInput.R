# analysis of PancreaticCancer dataset using iGC using Oncodrive Input
library(iGC)
rm(list = ls())

# load data
load("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset_OncodriveCISInput.Rdata")

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

# Computing gain/loss sample proportion ...
# 
# Computing CNA gain driven gene records ... 
# 
# |================================| 100%
# Computing CNA loss driven gene records ... 
# 
# |================================| 100%
# > toc()
# iGC: 15.142 sec elapsed


# saving results
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/")
write.table(cna_driven_genes$gain_driven, file = "iGC_PAAD_GainDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$loss_driven, file = "iGC_PAAD_LossDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$both, file = "iGC_PAAD_BothDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
