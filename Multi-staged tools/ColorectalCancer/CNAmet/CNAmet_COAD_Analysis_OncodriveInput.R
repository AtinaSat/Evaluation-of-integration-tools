# Analysis of ColonCancer dataset using CNAmet using OncodriveCIS input
library(CNAmet)
library(tictoc)
rm(list=ls())

# load data
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset_OncodriveCISInput.Rdata")

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


# Amplification: 1725.332 sec elapsed
# Deletion: 1680.907 sec elapsed
# Total: 3406.242 sec elapsed

ampResults <- data.frame(ampResults)
delResults <- data.frame(delResults)
ampResults$Gene <- rownames(ampResults)
delResults$Gene <- rownames(delResults)

# saving results
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/CNAmet/CNAmet_UsingOncodriveCISInput/")
write.table(ampResults, file = "CNAmet_COAD_GainDriven_Genes_OncodriveInput.tsv", row.names = F, sep = "\t", quote = F)
write.table(delResults, file = "CNAmet_COAD_LossDriven_Genes_OncodriveInput.tsv", row.names = F, sep = "\t", quote = F)




