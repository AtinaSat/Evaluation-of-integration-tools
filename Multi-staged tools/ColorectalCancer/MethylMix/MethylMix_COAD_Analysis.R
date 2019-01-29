# Analysis of Colorectal cancer dataset using MethylMix
# Since there are no normal samples, we just want the trancriptionally predictive genes. These can be obtained using 
# MethylMix_ModelGeneExpression
library(MethylMix)
library(doParallel)
library(tictoc)
rm(list = ls())

# load data
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset.Rdata")
rm(coad.cnv)

METCoad <- as.matrix(coad.me)
GECoad <- as.matrix(coad.ge)

# running methylMix in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
tic("MethylMix")
MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METCoad, GECoad, CovariateData = NULL)
toc()
stopCluster(cl)

# > MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METCoad, GECoad, CovariateData = NULL)
# Found 276 samples with both methylation and expression data.
# Correlating methylation data with gene expression...
# 
# Found 2156 transcriptionally predictive genes.
# > toc()
# MethylMix: 17.785 sec elapsed

# saving the results
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/MethylMix/")
write.table(MethylMixFunctionalGenes, "MethylMix_CompCOAD_3000Genes.txt", quote = F, row.names = F)
