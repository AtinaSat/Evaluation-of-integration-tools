## analysis of PancreaticCancer OncodriveIGC Input dataset using PLRS
library(plrs)

# load data
load("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset_OncodriveCISInput.Rdata")

# Reading GISTIC 2.0 copy number data
cnv_gistic <- read.table(gzfile("/home/anita/Integrated analysis in R/All_Cancers/PAAD_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"), 
                         header = T, sep = "\t")
rownames(cnv_gistic) <- cnv_gistic[,1]
cnv_gistic[,1] <- NULL
a <- intersect(rownames(paad.ge), rownames(cnv_gistic))
b <- intersect(colnames(paad.ge), colnames(cnv_gistic))
paad.cnv_gistic <- cnv_gistic[a,b]

# Substituting -2 to -1 in thresholded copy number data
paad.cnv[paad.cnv==-2] <- -1

paad.ge <- paad.ge[,order(colnames(paad.ge))]
paad.cnv <- paad.cnv[,order(colnames(paad.cnv))]
paad.cnv_gistic <- paad.cnv_gistic[,order(colnames(paad.cnv_gistic))]

all(colnames(paad.ge) == colnames(paad.cnv))
all(colnames(paad.ge) == colnames(paad.cnv_gistic))

paad.ge <- as.matrix(paad.ge)
paad.cnv <- as.matrix(paad.cnv)
paad.cnv_gistic <- as.matrix(paad.cnv_gistic)

# running plrs
library(tictoc)
tic("plrs")
resNoSel <- plrs.series(expr=paad.ge, cghseg=paad.cnv_gistic, cghcall=paad.cnv, control.select=NULL, control.model=list(min.obs=3))
toc()

# In progress... 
# 
# 10% done (1299 genes),	  time elapsed = 0:01:03
# 20% done (2596 genes),	  time elapsed = 0:02:05
# 30% done (3894 genes),	  time elapsed = 0:03:06
# 40% done (5191 genes),	  time elapsed = 0:04:11
# 50% done (6489 genes),	  time elapsed = 0:05:33
# 60% done (7787 genes),	  time elapsed = 0:07:18
# 70% done (9084 genes),	  time elapsed = 0:08:57
# 80% done (10382 genes),	  time elapsed = 0:10:10
# 90% done (11679 genes),	  time elapsed = 0:11:25
# 100% done (12977 genes),	  time elapsed = 0:12:40
# 
# > toc()
# plrs: 759.963 sec elapsed

summary(resNoSel)

# Results of test for each gene
head(resNoSel@test)
results <- data.frame(resNoSel@test)
results$Gene <- rownames(paad.ge)

# saving results
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/")
write.table(results, file = "PLRS_PAAD_Results_OncodriveInput.tsv", row.names = T, sep = "\t", quote = F)
