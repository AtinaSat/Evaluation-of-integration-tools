# creating simulated dataset as in the paper 
# 'Multi-omics integration - a comparison of unsupervised clustering methodologies'
# simulating case D no noise small dataset
# Group3 of Data type 1 and Group2 of Data type 2 and Group1 of Data type 3 not distinguishable
## Read breast cancer dataset 
# Paper: Comprehensive molecular portraits of human breast tumors

rm(list=ls())
gene_exp <- read.table("/home/anita/Benchmarking/molecular_portraits_breast/BRCA.exp.348.med.txt", header = T, sep = "\t", quote = "")
mirna_exp <- read.table("/home/anita/Benchmarking/molecular_portraits_breast/BRCA.348.mimat.txt", header = T, sep =",", quote = "")
me_exp <- read.table("/home/anita/Benchmarking/molecular_portraits_breast/BRCA.Methylation.574probes.802.txt", header = T, sep = "\t", quote = "")
metadata <- read.table("/home/anita/Benchmarking/molecular_portraits_breast/BRCA.sample.metedata.txt", header = T, sep = "\t", quote = "")

rownames(gene_exp) <- gene_exp[,1]
rownames(mirna_exp) <- mirna_exp[,1]
rownames(metadata) <- metadata[,1]

gene_exp[,1] <- NULL
mirna_exp[,1] <- NULL

colnames(gene_exp) <- substr(colnames(gene_exp),1,16)
colnames(mirna_exp) <- substr(colnames(mirna_exp),1,16)

gene_exp <- t(gene_exp)
mirna_exp <- t(mirna_exp)
me_exp <- t(me_exp)
head(colnames(gene_exp))

a <- intersect(rownames(gene_exp), rownames(mirna_exp))
b <- intersect(a, rownames(me_exp))
c <- intersect(b, rownames(metadata))

mrna <- gene_exp[c,]
mirna <- mirna_exp[c,]
meth <- me_exp[c,]
meta <- metadata[c,]

all(rownames(meth)==rownames(mrna))
all(rownames(meth)==rownames(mirna))
all(rownames(meta)==rownames(meth))

meta <- meta[,c(1,4)]

###############################
## for mixed signal
set.seed(999)
components <- sample(1:2,prob = c(0.5,0.5), size = 20, replace = T)

## for pca
source("/home/anita/utilities/pca_cor/PC-corr_net-master/R version/PC_corr_v2.R")
feat_names_ge <- c(paste0("Gene",as.character(1:500)))
feat_names <- c(paste0("Gene", as.character(1:500)))
sample_labels <- c(rep("Group1",20), rep("Group2",20), rep("Group3",20))
sample_names <- c(paste0("S",as.character(1:60)))

## simulation of gene expression data
geVar <- apply(mrna, 2, var)
Top1000 <- sort(geVar, decreasing = T)[1:1000]
mrna1000 <- mrna[, attr(Top1000, "names")]

geMeans <- apply(mrna1000, 2, mean)
geSds <- apply(mrna1000, 2, sd)

set.seed(999)
ind <- sample(1:1000, 500)
geSimulated <- matrix(,nrow = 60, ncol = 500)

for (i in 1:500){
  gmean <- geMeans[ind[i]]
  gsd <- geSds[ind[i]]
  # group1
  geSimulated[1:20,i] <- rnorm(n = 20, mean = gmean, sd = gsd) + rnorm(n=20, mean =0, sd = 0.4)
  # group2
  geSimulated[21:40,i] <- rnorm(n = 20, mean = gmean-(gmean/2), sd = gsd) + rnorm(n=20, mean =0, sd = 0.4)
  # group3
  mus = c(gmean, gmean-gmean/2)
  sds = c(gsd,gsd)
  geSimulated[41:60,i] <- rnorm(n = 20, mean = mus[components], sd = sds[components]) + rnorm(n=20, mean =0, sd = 0.4)
}

#PC_corr_v2(geSimulated,sample_labels,feat_names, sample_names)
# c("Group1","Group2","Group3")

colnames(geSimulated) <- feat_names_ge
rownames(geSimulated) <- sample_names

##################################################
## simulation of mirna expression data
miMeans <- apply(mirna, 2, mean)
miSds <- apply(mirna, 2, sd)

set.seed(444)
ind <- sample(1:dim(mirna)[2], 500)

mirnaSimulated <- matrix(,nrow = 60, ncol = 500)

for (i in 1:500){
  mirnamean <- miMeans[ind[i]]
  mirnasd <- miSds[ind[i]]
  # group1
  mirnaSimulated[1:20,i] <- rnorm(n = 20, mean = mirnamean, sd = mirnasd) + rnorm(n=20, mean =0, sd = 0.4)
  # group2
  mus = c(mirnamean, mirnamean+mirnamean/2)
  sds = c(mirnasd, mirnasd+(mirnasd/10))
  mirnaSimulated[21:40,i] <- rnorm(n = 20, mean = mus[components], sd = sds[components]) + rnorm(n=20, mean =0, sd = 0.4)
  # group3
  mirnaSimulated[41:60,i] <- rnorm(n = 20, mean = mirnamean+(mirnamean/2), sd = mirnasd+(mirnasd/10)) + rnorm(n=20, mean =0, sd = 0.4)
}

#PC_corr_v2(mirnaSimulated,sample_labels,feat_names, sample_names)

colnames(mirnaSimulated) <- feat_names
rownames(mirnaSimulated) <- sample_names

#################################################
## simulation of meth expression data
meMeans <- apply(meth, 2, mean)
meSds <- apply(meth, 2, sd)

set.seed(555)
ind <- sample(1:dim(meth)[2], 500, replace = T)
methSimulated <- matrix(,nrow = 60, ncol = 500)

for (i in 1:500){
  methmean <- meMeans[ind[i]]
  methsd <- meSds[ind[i]]
  # group1
  mus = c(methmean-methmean/2, methmean+methmean/2)
  sds = c(methsd, methsd+methsd/10)
  methSimulated[1:20,i] <- rnorm(n = 20, mean = mus[components], sd = sds[components]) + rnorm(n=20, mean =0, sd = 0.4)
  # group2
  methSimulated[21:40,i] <- rnorm(n = 20, mean = methmean-(methmean/2), sd = methsd) + rnorm(n=20, mean =0, sd = 0.4)
  # group3
  methSimulated[41:60,i] <- rnorm(n = 20, mean = methmean+(methmean/2), sd = methsd+(methsd/10)) + rnorm(n=20, mean =0, sd = 0.4)
}

#PC_corr_v2(methSimulated,sample_labels,feat_names, sample_names)

colnames(methSimulated) <- feat_names
rownames(methSimulated) <- sample_names

setwd("/home/anita/Benchmarking/simulated_data/")
save(geSimulated, mirnaSimulated, methSimulated, sample_labels, file= "SimulatedCaseD_NN.RData")
