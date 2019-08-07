# analysis of ColonCancer dataset using iGC using Oncodrive Input
library(iGC)
rm(list = ls())

# load data
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset_OncodriveCISInput.Rdata")

coad.cnv[coad.cnv==2] <- 1
coad.cnv[coad.cnv==-2] <- -1
coad.cnv <- cbind(rownames(coad.cnv), coad.cnv)
colnames(coad.cnv)[1] <- "GENE"
coad.cnv[,1] <- as.character(coad.cnv[,1])

coad.ge <- cbind(rownames(coad.ge), coad.ge)
colnames(coad.ge)[1] <- "GENE"
coad.ge[,1] <- as.character(coad.ge[,1])

# Change the data types to comply with the tool requirement
library(data.table)
setDT(coad.cnv)
setDT(coad.ge)
class(coad.cnv)
class(coad.ge)

# running iGC
library(tictoc)
tic("iGC")
cna_driven_genes <- find_cna_driven_gene(coad.cnv, coad.ge)
toc()

# Computing gain/loss sample proportion ...
# 
# Computing CNA gain driven gene records ... 
# 
# |================================| 100%
# Computing CNA loss driven gene records ... 
# 
# |================================| 100%
# 
# iGC: 21.213 sec elapsed

# saving results
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/")
write.table(cna_driven_genes$gain_driven, file = "iGC_COAD_GainDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$loss_driven, file = "iGC_COAD_LossDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$both, file = "iGC_COAD_BothDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
