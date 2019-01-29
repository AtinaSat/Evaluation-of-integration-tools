## analysis of pancreatic cancer dataset using PLRS
library(plrs)

# load data
load("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/PancreaticCancerRawDataset.Rdata")
rm(paad.me)

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

# > resNoSel <- plrs.series(expr=paad.ge, cghseg=paad.cnv_gistic, cghcall=paad.cnv, control.select=NULL, control.model=list(min.obs=3))
# In progress... 
# 
# 10% done (1615 genes),	  time elapsed = 0:01:19
# 20% done (3230 genes),	  time elapsed = 0:02:42
# 30% done (4844 genes),	  time elapsed = 0:04:05
# 40% done (6458 genes),	  time elapsed = 0:05:20
# 50% done (8072 genes),	  time elapsed = 0:06:40
# 60% done (9687 genes),	  time elapsed = 0:07:55
# 70% done (11301 genes),	  time elapsed = 0:09:10
# 80% done (12915 genes),	  time elapsed = 0:10:26
# 90% done (14530 genes),	  time elapsed = 0:11:40
# 100% done (16144 genes),	  time elapsed = 0:13:03
# 
# > toc()
# plrs: 782.999 sec elapsed

summary(resNoSel)

# Results of test for each gene
head(resNoSel@test)
results <- data.frame(resNoSel@test)
results$Gene <- rownames(paad.ge)

# saving results
setwd("/home/anita/Benchmarking/two_omics/PancreaticCancerCompleteDataAnalysis/plrs")
write.table(results, file = "PLRS_PAAD_Results.tsv", row.names = T, sep = "\t", quote = F)
