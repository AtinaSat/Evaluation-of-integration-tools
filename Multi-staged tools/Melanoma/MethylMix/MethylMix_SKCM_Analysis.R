# MethylMix analysis of Melanoma dataset
# Since there are no normal samples, we just want the trancriptionally predictive genes. These can be obtained using 
# MethylMix_ModelGeneExpression

library(MethylMix)
library(doParallel)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaRawDataset.Rdata")
source("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MethylMix/MethylMix_FDRCalculation_Function.R")
rm(skcm.cnv)

METSkcm <- as.matrix(skcm.me)
GESkcm <- as.matrix(skcm.ge)

# running methylmix in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
tic("MethylMix")
MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METSkcm, GESkcm, CovariateData = NULL)
toc()
MethylMixFunctionalGenes_FDR <- MethylMix_ModelGeneExpression_FDR(METcancer = METSkcm, GEcancer = GESkcm, CovariateData = NULL)
stopCluster(cl)

# Found 366 samples with both methylation and expression data.
# Correlating methylation data with gene expression...
# 
# Found 2492 transcriptionally predictive genes.
# > toc()
# MethylMix: 15.47 sec elapsed

fdrs <- data.frame(MethylMixFunctionalGenes_FDR)
fdrs[,2] <- as.numeric(as.character(fdrs[,2]))
str(fdrs)

fdrs2 <- fdrs[fdrs[,2] < 0.05,] # All the genes satisfy fdr cut-off of 0.05
genes <- as.character(MethylMixFunctionalGenes)

# saving results
setwd("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MethylMix/")
write.table(genes, "MethylMix_SKCM_3000Genes.txt", quote = F, row.names = F)
