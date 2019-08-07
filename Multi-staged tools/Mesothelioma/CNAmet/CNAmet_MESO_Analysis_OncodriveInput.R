# Analysis of Mesothelioma dataset using CNAmet using OncodriveCIS input
library(CNAmet)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset_OncodriveCISInput.Rdata")

### reading copy number values and creating two binary matrices for amplification and deletion calls
cn.amp <- meso.cnv
cn.amp[cn.amp>1] <- 1
cn.amp[cn.amp<0] <- 0

cn.del <- meso.cnv
cn.del[cn.del>0] <- 0
cn.del[cn.del<0] <- 1
cn.del[1:5,1:5]

cn.amp <- as.matrix(cn.amp, rownames.force = T)
cn.del <- as.matrix(cn.del, rownames.force = T)
ge.exp <- as.matrix(meso.ge, rownames.force = T)

# running CNAmet
tic("Total")
tic("Amplification")
ampResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = cn.amp, methylMatrix = NULL, perms = 1000, na.limit = 0.1, gainData = TRUE,
                     favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
toc()
tic("Deletion")
delResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = cn.del, methylMatrix = NULL, perms = 1000, na.limit = 0.1, gainData = FALSE,
                     favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
toc()
toc()


ampResults <- data.frame(ampResults)
delResults <- data.frame(delResults)
ampResults$Gene <- rownames(ampResults)
delResults$Gene <- rownames(delResults)

# saving results
setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/")
write.table(ampResults, file = "CNAmet_MESO_GainDriven_Genes_OncodriveInput.tsv", row.names = F, sep = "\t", quote = F)
write.table(delResults, file = "CNAmet_MESO_LossDriven_Genes_OncodriveInput.tsv", row.names = F, sep = "\t", quote = F)




