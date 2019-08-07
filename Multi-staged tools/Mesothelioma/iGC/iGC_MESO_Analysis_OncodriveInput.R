# analysis of Mesothelioma dataset using iGC using Oncodrive Input
library(iGC)
rm(list = ls())

# load data
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset_OncodriveCISInput.Rdata")

meso.cnv[meso.cnv==2] <- 1
meso.cnv[meso.cnv==-2] <- -1
meso.cnv <- cbind(rownames(meso.cnv), meso.cnv)
colnames(meso.cnv)[1] <- "GENE"
meso.cnv[,1] <- as.character(meso.cnv[,1])

meso.ge <- cbind(rownames(meso.ge), meso.ge)
colnames(meso.ge)[1] <- "GENE"
meso.ge[,1] <- as.character(meso.ge[,1])

# Change the data types to comply with the tool requirement
library(data.table)
setDT(meso.cnv)
setDT(meso.ge)
class(meso.cnv)
class(meso.ge)

# running iGC
library(tictoc)
tic("iGC")
cna_driven_genes <- find_cna_driven_gene(meso.cnv, meso.ge)
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
# iGC: 10.023 sec elapsed

# saving results
setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/iGC/iGC_UsingOncodriveCISInput/")
write.table(cna_driven_genes$gain_driven, file = "iGC_MESO_GainDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$loss_driven, file = "iGC_MESO_LossDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$both, file = "iGC_MESO_BothDriven_genes_OncodriveInput.tsv", sep = "\t", quote = F,row.names = F )
