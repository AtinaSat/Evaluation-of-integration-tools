# Analysis of colorectal cancer dataset with CNAmet
library(CNAmet)
library(mltools)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset.Rdata")
rm(coad.me)

### reading copy number values and creating two binary matrices for amplification and deletion calls
cn.amp <- coad.cnv
cn.amp[cn.amp>1] <- 1
cn.amp[cn.amp<0] <- 0

cn.del <- coad.cnv
cn.del[cn.del>0] <- 0
cn.del[cn.del<0] <- 1
cn.del[1:5,1:5]

cn.amp <- as.matrix(cn.amp, rownames.force = T)
cn.del <- as.matrix(cn.del, rownames.force = T)
ge.exp <- as.matrix(coad.ge, rownames.force = T)

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

# > ampResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = cn.amp, methylMatrix = NULL, perms = 1000, na.limit = 0.1, gainData = TRUE,
#                        +                      favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
# > toc()
# Amplification: 2097.784 sec elapsed
# > tic("Deletion")
# > delResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = cn.del, methylMatrix = NULL, perms = 1000, na.limit = 0.1, gainData = FALSE,
#                        +                      favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
# > toc()
# Deletion: 2115.862 sec elapsed

ampResults <- data.frame(ampResults)
delResults <- data.frame(delResults)
ampResults$Gene <- rownames(ampResults)
delResults$Gene <- rownames(delResults)

# saving results
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/")
write.table(ampResults, file = "CNAmet_CompCOAD_GainDriven_Genes.tsv", row.names = F, sep = "\t", quote = F)
write.table(delResults, file = "CNAmet_CompCOAD_LossDriven_Genes.tsv", row.names = F, sep = "\t", quote = F)

