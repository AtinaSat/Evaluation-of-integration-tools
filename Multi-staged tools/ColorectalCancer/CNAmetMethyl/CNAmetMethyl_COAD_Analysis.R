# Analysis of colorectal cancer dataset using CNAmet - Methylation and gene expression integration
library(CNAmet)
library(tictoc)
rm(list=ls())

load("../ColorectalCancerRawDataset.Rdata")

## binning the methylation values into 3 bins (hypo, normal, and hyper)
me.exp.binned <- coad.me
hist(me.exp.binned[,1], xlim = c(0,1))
me.exp.binned[1:5,1:5]
me.exp.binned[me.exp.binned>0.8] <- 1
me.exp.binned[me.exp.binned<0.2] <- -1
me.exp.binned[me.exp.binned!= -1 & me.exp.binned!= 1] <- 0

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
ge.exp <- as.matrix(coad.ge, rownames.force = T)

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

# tic("Total")
# > tic("Hypomethylation")
# > hypoResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = NULL, methylMatrix = me.hypo, perms = 1000, na.limit = 0.1, gainData = TRUE,
#                         +                      favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
# > toc()
# Hypomethylation: 1405.844 sec elapsed
# > tic("Hypermethylation")
# > hyperResults <- CNAmet(exprMatrix = ge.exp, cghMatrix = NULL, methylMatrix = me.hyper, perms = 1000, na.limit = 0.1, gainData = FALSE,
#                        +                      favorSynergetic = TRUE, strictChecks = FALSE, strictLim = 0.05)
# > toc()
# Hypermethylation: 1279.988 sec elapsed
# > toc()
# Total: 2685.834 sec elapsed

hypoResults <- data.frame(hypoResults)
hyperResults <- data.frame(hyperResults)
hypoResults$Gene <- rownames(hypoResults)
hyperResults$Gene <- rownames(hyperResults)

# saving results
write.table(hypoResults, file = "CNAmetMethyl_COAD_HypoDriven_Genes.tsv", row.names = F, sep = "\t", quote = F)
write.table(hyperResults, file = "CNAmetMethyl_COAD_HyperDriven_Genes.tsv", row.names = F, sep = "\t", quote = F)

