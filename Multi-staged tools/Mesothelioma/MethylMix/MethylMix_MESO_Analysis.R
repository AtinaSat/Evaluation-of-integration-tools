# MethylMix analysis of Mesothelioma dataset
# Since there are no normal samples, we just want the trancriptionally predictive genes. These can be obtained using 
# MethylMix_ModelGeneExpression

library(MethylMix)
library(doParallel)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset.Rdata")
source("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MethylMix/MethylMix_FDRCalculation_Function.R")
rm(meso.cnv)

METMeso <- as.matrix(meso.me)
GEMeso <- as.matrix(meso.ge)

# running methylmix in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
tic("MethylMix")
MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METcancer = METMeso, GEcancer = GEMeso, CovariateData = NULL)
toc()
MethylMixFunctionalGenes_FDR <- MethylMix_ModelGeneExpression_FDR(METcancer = METMeso, GEcancer = GEMeso, CovariateData = NULL)
stopCluster(cl)

# Found 87 samples with both methylation and expression data.
# Correlating methylation data with gene expression...
# 
# Found 2346 transcriptionally predictive genes.
# > toc()
# MethylMix: 15.929 sec elapsed

fdrs <- data.frame(MethylMixFunctionalGenes_FDR)
fdrs[,2] <- as.numeric(as.character(fdrs[,2]))
str(fdrs)

fdrs2 <- fdrs[fdrs[,2] < 0.05,] # All the genes satisfy fdr cut-off of 0.05

genes <- as.character(MethylMixFunctionalGenes)

# saving results
setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MethylMix/")
write.table(genes, "MethylMix_MESO_3000Genes.txt", quote = F, row.names = F)
