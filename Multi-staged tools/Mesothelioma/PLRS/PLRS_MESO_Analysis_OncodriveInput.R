## analysis of Mesothelioma cancer dataset using PLRS
library(plrs)

# load data
load("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaRawDataset_OncodriveCISInput.Rdata")

# Reading GISTIC 2.0 copy number data
cnv_gistic <- read.table(gzfile("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/MesotheliomaData/MESO_Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"), 
                         header = T, sep = "\t")
rownames(cnv_gistic) <- cnv_gistic[,1]
cnv_gistic[,1] <- NULL
a <- intersect(rownames(meso.ge), rownames(cnv_gistic))
b <- intersect(colnames(meso.ge), colnames(cnv_gistic))
meso.cnv_gistic <- cnv_gistic[a,b]

# Substituting -2 to -1 in thresholded copy number data
meso.cnv[meso.cnv==-2] <- -1

meso.ge <- meso.ge[,order(colnames(meso.ge))]
meso.cnv <- meso.cnv[,order(colnames(meso.cnv))]
meso.cnv_gistic <- meso.cnv_gistic[,order(colnames(meso.cnv_gistic))]

all(colnames(meso.ge) == colnames(meso.cnv))
all(colnames(meso.ge) == colnames(meso.cnv_gistic))

meso.ge <- as.matrix(meso.ge)
meso.cnv <- as.matrix(meso.cnv)
meso.cnv_gistic <- as.matrix(meso.cnv_gistic)

# running plrs
library(tictoc)
tic("plrs")
resNoSel <- plrs.series(expr=meso.ge, cghseg=meso.cnv_gistic, cghcall=meso.cnv, control.select=NULL, control.model=list(min.obs=3))
toc()

# In progress... 
# 
# 10% done (838 genes),	  time elapsed = 0:00:18
# 20% done (1675 genes),	  time elapsed = 0:00:36
# 30% done (2512 genes),	  time elapsed = 0:00:54
# 40% done (3349 genes),	  time elapsed = 0:01:11
# 50% done (4186 genes),	  time elapsed = 0:01:29
# 60% done (5024 genes),	  time elapsed = 0:01:47
# 70% done (5861 genes),	  time elapsed = 0:02:05
# 80% done (6698 genes),	  time elapsed = 0:02:22
# 90% done (7535 genes),	  time elapsed = 0:02:38
# 100% done (8372 genes),	  time elapsed = 0:02:56
# 
# Warning message:
#   In plrs(expr = expr[x, ], cghseg = cghseg[x, ], cghcall = cghcall[x,  :
#                                                                       No gains: observations labelled amplifications are changed to gains
# > toc()
# plrs: 176.263 sec elapsed

summary(resNoSel)

# Results of test for each gene
head(resNoSel@test)
results <- data.frame(resNoSel@test)
results$Gene <- rownames(meso.ge)

# saving results
setwd("/home/anita/Benchmarking/two_omics/MesotheliomaCompleteDataAnalysis/plrs/PLRS_UsingOncodriveCISInput/")
write.table(results, file = "PLRS_MESO_Results_OncodriveInput.tsv", row.names = T, sep = "\t", quote = F)
