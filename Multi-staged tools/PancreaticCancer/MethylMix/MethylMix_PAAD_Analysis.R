# MethylMix analysis of pancreatic cancer dataset
# Since there are no normal samples, we just want the trancriptionally predictive genes. These can be obtained using 
# MethylMix_ModelGeneExpression

library(MethylMix)
library(doParallel)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset.Rdata")
source("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MethylMix/MethylMix_FDRCalculation_Function.R")
rm(meso.cnv)
rm(paad.cnv)

METPaad <- as.matrix(paad.me)
GEPaad <- as.matrix(paad.ge)

# running methylmix in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
tic("MethylMix")
MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METPaad, GEPaad, CovariateData = NULL)
toc()
MethylMixFunctionalGenes_FDR <- MethylMix_ModelGeneExpression_FDR(METcancer = METPaad, GEcancer = GEPaad, CovariateData = NULL)
stopCluster(cl)

# > MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METPaad, GEPaad, CovariateData = NULL)
# Found 177 samples with both methylation and expression data.
# Correlating methylation data with gene expression...
# 
# Found 3025 transcriptionally predictive genes.
# > toc()
# MethylMix: 16.389 sec elapsed

fdrs <- data.frame(MethylMixFunctionalGenes_FDR)
fdrs[,2] <- as.numeric(as.character(fdrs[,2]))
str(fdrs)

fdrs2 <- fdrs[fdrs[,2] < 0.05,] # All the genes satisfy fdr cut-off of 0.05
# sorting based on fdr (smallest to largest)
head(fdrs2)
fdrs2 <- fdrs2[order(fdrs2$V2, decreasing = F),]
head(fdrs2)

# selecting top 3000
resdf <- fdrs2[1:3000,]
genes <- as.character(resdf$FunctionalGenes)

# saving results
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/MethylMix/")
write.table(genes, "MethylMix_PAAD_3000Genes.txt", quote = F, row.names = F)
