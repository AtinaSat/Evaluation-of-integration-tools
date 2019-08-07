## analysis of Melanoma cancer dataset using PLRS
library(plrs)

# load data
load("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaRawDataset_OncodriveCISInput.Rdata")

# Reading GISTIC 2.0 copy number data
cnv_gistic <- read.table(gzfile("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/MelanomaData/SKCM_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"), 
                         header = T, sep = "\t")
rownames(cnv_gistic) <- cnv_gistic[,1]
cnv_gistic[,1] <- NULL
a <- intersect(rownames(skcm.ge), rownames(cnv_gistic))
b <- intersect(colnames(skcm.ge), colnames(cnv_gistic))
skcm.cnv_gistic <- cnv_gistic[a,b]

# Substituting -2 to -1 in thresholded copy number data
skcm.cnv[skcm.cnv==-2] <- -1

skcm.ge <- skcm.ge[,order(colnames(skcm.ge))]
skcm.cnv <- skcm.cnv[,order(colnames(skcm.cnv))]
skcm.cnv_gistic <- skcm.cnv_gistic[,order(colnames(skcm.cnv_gistic))]

all(colnames(skcm.ge) == colnames(skcm.cnv))
all(colnames(skcm.ge) == colnames(skcm.cnv_gistic))

skcm.ge <- as.matrix(skcm.ge)
skcm.cnv <- as.matrix(skcm.cnv)
skcm.cnv_gistic <- as.matrix(skcm.cnv_gistic)

# running plrs
library(tictoc)
tic("plrs")
resNoSel <- plrs.series(expr=skcm.ge, cghseg=skcm.cnv_gistic, cghcall=skcm.cnv, control.select=NULL, control.model=list(min.obs=3))
toc()
# In progress... 
# 
# 10% done (1475 genes),	  time elapsed = 0:01:50
# 20% done (2950 genes),	  time elapsed = 0:03:39
# 30% done (4424 genes),	  time elapsed = 0:05:34
# 40% done (5899 genes),	  time elapsed = 0:07:29
# 50% done (7373 genes),	  time elapsed = 0:09:22
# 60% done (8847 genes),	  time elapsed = 0:11:12
# 70% done (10322 genes),	  time elapsed = 0:13:04
# 80% done (11796 genes),	  time elapsed = 0:14:53
# 90% done (13271 genes),	  time elapsed = 0:16:42
# 100% done (14745 genes),	  time elapsed = 0:18:31
# 
# > toc()
# plrs: 1111.408 sec elapsed

summary(resNoSel)

# Results of test for each gene
head(resNoSel@test)
results <- data.frame(resNoSel@test)
results$Gene <- rownames(skcm.ge)

# saving results
setwd("/home/anita/Benchmarking/two_omics/MelanomaCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/")
write.table(results, file = "PLRS_SKCM_Results_OncodriveInput.tsv", row.names = T, sep = "\t", quote = F)
