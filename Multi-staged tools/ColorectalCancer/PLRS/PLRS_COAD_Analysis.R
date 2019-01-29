## Analysis of colorectal cancer dataset using PLRS
library(plrs)
rm(list = ls())

#load data
load("../ColonrectalCancerRawDataset.Rdata")
rm(coad.me)

# Reading GISTIC 2.0 copy number data
cnv_gistic <- read.table(gzfile("COAD_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"), 
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

# > resNoSel <- plrs.series(expr=coad.ge, cghseg=coad.cnv_gistic, cghcall=coad.cnv, control.select=NULL, control.model=list(min.obs=3))
# In progress... 
# 
# 10% done (1594 genes),	  time elapsed = 0:00:45
# 20% done (3188 genes),	  time elapsed = 0:01:33
# 30% done (4781 genes),	  time elapsed = 0:02:28
# 40% done (6375 genes),	  time elapsed = 0:03:15
# 50% done (7968 genes),	  time elapsed = 0:04:06
# 60% done (9561 genes),	  time elapsed = 0:04:60
# 70% done (11155 genes),	  time elapsed = 0:05:49
# 80% done (12748 genes),	  time elapsed = 0:06:38
# 90% done (14342 genes),	  time elapsed = 0:07:28
# 100% done (15935 genes),	  time elapsed = 0:08:16
# 
# > toc()
# plrs: 496.392 sec elapsed

summary(resNoSel)

# Results of test for each gene
head(resNoSel@test)
results <- data.frame(resNoSel@test)
results$Gene <- rownames(coad.ge)

# saving results
write.table(results, file = "PLRS_COAD_Results.tsv", row.names = T, sep = "\t", quote = F)



