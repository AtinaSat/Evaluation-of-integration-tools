#' ---
#' title: "GBM data"
#' author: "Anita"
#' date: "27 June 2019"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ------------------------------------------------------------------------
library(SNFtool)
library(tictoc)
gbm.dat.me = readRDS("GBM_Methylation_Top25Filtered.rd")
gbm.dat.rna = readRDS("GBM_RNA_Top25Filtered.rd")
gbm.dat.protein = readRDS("GBM_Protein_Top25Filtered.rd")
data.labels = readRDS("GBM_DataLabels.rd")

# Since the are normalized we do not perform further normalization
tic("GBM 4Labels")
# Calculating distance matrices
dist1 <- as.matrix(dist(gbm.dat.rna))
dist2 <- as.matrix(dist(gbm.dat.protein))
dist3 <- as.matrix(dist(gbm.dat.me))

# Calculating affinity matrics
W1 <- affinityMatrix(dist1)
W2 <- affinityMatrix(dist2)
W3 <- affinityMatrix(dist3)

# Combining the clusters
W = SNF(list(W1,W2, W3))

# Spectral clustering
clusters_label <- data.frame(spectralClustering(W,K = 4))

toc()

# GBM 4Labels: 0.616 sec elapsed


## Clustering results
clusters_label$Samples <- data.labels$Label
clusters_label$Samples <- as.character(clusters_label$Samples)
table(clusters_label$spectralClustering.W..K...4.)

clusters_label$Truth1 <- clusters_label$Samples
confusion.mat = table(true = clusters_label$Truth1, Pred = clusters_label$spectralClustering.W..K...4.)
confusion.mat
# Pred
# true           1  2  3  4
# Classical   23  0  5  0
# Mesenchymal  2  0 17  5
# Neural       6  3  8  0
# Proneural    2 23  1  5


## CLASSICAL: CLUSTER 1; MESENCHYMAL: CLUSTER 3; NEURAL: CLUSTER 4; PRONEURAL: CLUSTER 2

clusters_label$Truth <- clusters_label$Samples
clusters_label$Truth <- gsub("Classical", "1", clusters_label$Truth)
clusters_label$Truth <- gsub("Mesenchymal", "3", clusters_label$Truth)
clusters_label$Truth <- gsub("Neural", "4", clusters_label$Truth)
clusters_label$Truth <- gsub("Proneural", "2", clusters_label$Truth)
head(clusters_label)
confusion.mat = table(true = clusters_label$Truth, Pred = clusters_label$spectralClustering.W..K...4.)
confusion.mat
# Pred
# true  1  2  3  4
# 1 23  0  5  0
# 2  2 23  1  5
# 3  2  0 17  5
# 4  6  3  8  0

#' 
#' Performance analysis
#' 
## ------------------------------------------------------------------------
SNF_Accuracy = sum(clusters_label$Truth == clusters_label$spectralClustering.W..K...4.)/length(clusters_label$spectralClustering.W..K...4.)
round(SNF_Accuracy,9)
# 0.63

# Classical cluster 1
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0

# PRONEURAL Cluster 2
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# MESENCHYMAL Cluster 3
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0

# NEURAL Cluster 4
Pg4 <- round(confusion.mat[4,4]/sum(confusion.mat[,4]),2) # Precision
Rg4 <- round(confusion.mat[4,4]/sum(confusion.mat[4,]),2) # Recall
Fg4 <- round(2*((Pg4*Rg4)/(Pg4+Rg4)),2) # F1-score
Fg4[is.nan(Fg4)] <- 0

AvgF <- round((Fg1+Fg2+Fg3+Fg4)/4,2)
# 66
#' 
#' 
## ------------------------------------------------------------------------
setwd("/home/anita/Benchmarking/GBM/SNF")
save.image(file = "GBM_SNFTool_Results.RData")

