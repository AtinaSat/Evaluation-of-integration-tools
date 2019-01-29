# analysis of colorectal cancer dataset using iGC
library(iGC)
rm(list = ls())

# load data
load("../ColorectalCancerRawDataset.Rdata")
rm(coad.me)

# Data processing to fit tool input
coad.cnv[coad.cnv==2] <- 1
coad.cnv[coad.cnv==-2] <- -1
coad.cnv <- cbind(rownames(coad.cnv), coad.cnv)
colnames(coad.cnv)[1] <- "GENE"
coad.cnv[,1] <- as.character(coad.cnv[,1])

coad.ge <- cbind(rownames(coad.ge), coad.ge)
colnames(coad.ge)[1] <- "GENE"
coad.ge[,1] <- as.character(coad.ge[,1])

# Change the data types to comply with the tool requirement
library(data.table)
setDT(coad.cnv)
setDT(coad.ge)
class(coad.cnv)
class(coad.ge)

# running iGC
library(tictoc)
tic("iGC")
cna_driven_genes <- find_cna_driven_gene(coad.cnv, coad.ge)
toc()

# > cna_driven_genes <- find_cna_driven_gene(coad.cnv, coad.ge)
# Computing gain/loss sample proportion ...
# 
# Computing CNA gain driven gene records ... 
# 
# |================================| 100%
# Computing CNA loss driven gene records ... 
# 
# |================================| 100%
# > toc()
# iGC: 23.844 sec elapsed

# saving results
write.table(cna_driven_genes$gain_driven, file = "iGC_COAD_GainDriven_genes.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$loss_driven, file = "iGC_COAD_LossDriven_genes.tsv", sep = "\t", quote = F,row.names = F )
write.table(cna_driven_genes$both, file = "iGC_COAD_BothDriven_genes.tsv", sep = "\t", quote = F,row.names = F )
