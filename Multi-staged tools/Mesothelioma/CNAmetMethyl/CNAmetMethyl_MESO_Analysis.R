# Analysis if Mesothelioma cancer dataset using CNAmet - methylation and gene expression integration
library(CNAmet)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset.Rdata")

## binning the methylation values into 3 bins (hypo, normal, and hyper)
me.exp.binned <- meso.me
hist(me.exp.binned[,1], xlim = c(0,1))
me.exp.binned[1:5,1:5]
me.exp.binned[me.exp.binned>0.8] <- 1 # hypermethylation
me.exp.binned[me.exp.binned<0.2] <- -1 # hypomethylation
me.exp.binned[me.exp.binned!= -1 & me.exp.binned!= 1] <- 0 # normal

table(me.exp.binned[,1])

# creating binary matrices for hypo and hypermethylation data
me.hyper <- me.exp.binned 
table(me.hyper[,1])
me.hyper[me.hyper<1] <-0 
me.hyper[me.hyper>0] <- 1

me.hypo <- me.exp.binned
table(me.hypo[,1])
me.hypo[me.hypo>-1] <- 0
me.hypo[me.hypo<0] <- 1

me.hyper <- as.matrix(me.hyper, rownames.force = T)
me.hypo <- as.matrix(me.hypo, rownames.force = T)
ge.exp <- as.matrix(meso.ge, rownames.force = T)

# running CNAmet
tic("Total")
tic("Hypomethylation")
hypoResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = NULL, methylMatrix = me.hypo, perms = 1000, na.limit = 0.1, gainData = TRUE,
                      favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
toc()
tic("Hypermethylation")
hyperResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = NULL, methylMatrix = me.hyper, perms = 1000, na.limit = 0.1, gainData = FALSE,
                       favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
toc()
toc()

# Hypomethylation: 1070.401 sec elapsed
# Hypermethylation: 1036.054 sec elapsed
# Total: 2106.457 sec elapsed

hypoResults <- data.frame(hypoResults)
hyperResults <- data.frame(hyperResults)
hypoResults$Gene <- rownames(hypoResults)
hyperResults$Gene <- rownames(hyperResults)

# saving results
setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/CNAmetMethyl/")
write.table(hypoResults, file = "CNAmetMethyl_MESO_HypoDriven_Genes.tsv", row.names = F, sep = "\t", quote = F)
write.table(hyperResults, file = "CNAmetMethyl_MESO_HyperDriven_Genes.tsv", row.names = F, sep = "\t", quote = F)

