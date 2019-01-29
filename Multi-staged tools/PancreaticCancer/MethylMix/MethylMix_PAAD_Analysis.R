# MethylMix analysis of pancreatic cancer dataset
# Since there are no normal samples, we just want the trancriptionally predictive genes. These can be obtained using 
# MethylMix_ModelGeneExpression

library(MethylMix)
library(doParallel)
library(tictoc)
rm(list=ls())

# load data
load("PancreaticCancerRawDataset.Rdata")
rm(paad.cnv)

METPaad <- as.matrix(paad.me)
GEPaad <- as.matrix(paad.ge)

# running methylmix in parallel
cl <- makeCluster(5)
registerDoParallel(cl)
tic("MethylMix")
MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METPaad, GEPaad, CovariateData = NULL)
toc()
stopCluster(cl)

# > MethylMixFunctionalGenes <- MethylMix_ModelGeneExpression(METPaad, GEPaad, CovariateData = NULL)
# Found 177 samples with both methylation and expression data.
# Correlating methylation data with gene expression...
# 
# Found 3025 transcriptionally predictive genes.
# > toc()
# MethylMix: 16.389 sec elapsed

# Since the output is just a list of genes, we rank the genes based on the pearson correlation
# between the methylation and gene expression of the gene and then select the top 3000

METPaad2 <- METPaad[MethylMixFunctionalGenes,]
GEPaad2 <- GEPaad[MethylMixFunctionalGenes,]

rhoVals <- sapply(1:nrow(METPaad2), function(i) cor(METPaad2[i,], GEPaad2[i,]))

resdf <- data.frame(MethylMixFunctionalGenes)
resdf$RhoVals <- rhoVals
resdf <- resdf[order(resdf$RhoVals),]
resdf <- resdf[1:3000,]

genes <- as.character(resdf$MethylMixFunctionalGenes)

# saving results
write.table(genes, "MethylMix_PAAD_3000Genes.txt", quote = F, row.names = F)
