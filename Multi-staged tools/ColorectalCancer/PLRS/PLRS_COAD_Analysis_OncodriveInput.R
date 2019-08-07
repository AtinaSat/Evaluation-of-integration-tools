## analysis of ColonCancer OncodriveCIS Input dataset using PLRS
library(plrs)

# load data
load("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/ColonCancerRawDataset_OncodriveCISInput.Rdata")

# Reading GISTIC 2.0 copy number data
cnv_gistic <- read.table(gzfile("/home/anita/Integrated analysis in R/All_Cancers/COAD_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"), 
                         header = T, sep = "\t")
rownames(cnv_gistic) <- cnv_gistic[,1]
cnv_gistic[,1] <- NULL
a <- intersect(rownames(coad.ge), rownames(cnv_gistic))
b <- intersect(colnames(coad.ge), colnames(cnv_gistic))
coad.cnv_gistic <- cnv_gistic[a,b]

# Substituting -2 to -1 in thresholded copy number data
coad.cnv[coad.cnv==-2] <- -1

coad.ge <- coad.ge[,order(colnames(coad.ge))]
coad.cnv <- coad.cnv[,order(colnames(coad.cnv))]
coad.cnv_gistic <- coad.cnv_gistic[,order(colnames(coad.cnv_gistic))]

all(colnames(coad.ge) == colnames(coad.cnv))
all(colnames(coad.ge) == colnames(coad.cnv_gistic))

coad.ge <- as.matrix(coad.ge)
coad.cnv <- as.matrix(coad.cnv)
coad.cnv_gistic <- as.matrix(coad.cnv_gistic)

# running plrs
library(tictoc)
tic("plrs")
resNoSel <- plrs.series(expr=coad.ge, cghseg=coad.cnv_gistic, cghcall=coad.cnv, control.select=NULL, control.model=list(min.obs=3))
toc()

# In progress... 
# 
# 10% done (1237 genes),	  time elapsed = 0:00:39
# 20% done (2472 genes),	  time elapsed = 0:01:16
# 30% done (3708 genes),	  time elapsed = 0:01:55
# 40% done (4944 genes),	  time elapsed = 0:02:36
# 50% done (6180 genes),	  time elapsed = 0:03:17
# 60% done (7415 genes),	  time elapsed = 0:04:02
# 70% done (8651 genes),	  time elapsed = 0:04:49
# 80% done (9887 genes),	  time elapsed = 0:05:31
# 90% done (11122 genes),	  time elapsed = 0:06:12
# 100% done (12358 genes),	  time elapsed = 0:06:46
# 
# > toc()
# plrs: 406.489 sec elapsed


summary(resNoSel)

# Results of test for each gene
head(resNoSel@test)
results <- data.frame(resNoSel@test)
results$Gene <- rownames(coad.ge)

# saving results
setwd("/home/anita/Benchmarking/two_omics/ColonCancerCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/")
write.table(results, file = "PLRS_COAD_Results_OncodriveInput.tsv", row.names = T, sep = "\t", quote = F)
