#' ---
#' title: "Hepatocellular Carcinoma  SNF - 3 labels"
#' author: "Anita"
#' date: "7 September 2018"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ------------------------------------------------------------------------
library(SNFtool)
library(tictoc)
load("/home/anita/Benchmarking/HepatocellularCarcinoma/data/HPC_filtered_data.Rdata")
mirna2 <- log(mirna+1)
mrna2 <- log(mrna+1)

#' 
#' performing SNF
#' 
## ------------------------------------------------------------------------

tic("SNF-HPC3Labels")

# Calculating distance matrices
dist1 <- as.matrix(dist(mrna2))
dist2 <- as.matrix(dist(mirna2))
dist3 <- as.matrix(dist(me))

# Calculating affinity matrics
W1 <- affinityMatrix(dist1)
W2 <- affinityMatrix(dist2)
W3 <- affinityMatrix(dist3)

# Combining the clusters
W = SNF(list(W1,W2, W3))

# Spectral clustering
clusters_label <- data.frame(spectralClustering(W,K = 3))

toc()

#' 
#' ## Clustering results
#' 
## ------------------------------------------------------------------------
clusters_label$Samples <- data.labels$Samples
clusters_label$Samples <- as.character(clusters_label$Samples)


clusters_label$Truth <- clusters_label$Samples
clusters_label$Truth <- gsub("N.{1,4}", "2", clusters_label$Truth)
clusters_label$Truth <- gsub("T.{1,4}", "1", clusters_label$Truth)
clusters_label$Truth <- gsub("P.{1,4}", "3", clusters_label$Truth)

confusion.mat = table(true = clusters_label$Truth, Pred = clusters_label$spectralClustering.W..K...3.)
confusion.mat

#' 
#' Performance analysis
#' 
## ------------------------------------------------------------------------
SNF_Accuracy = sum(clusters_label$Truth == clusters_label$spectralClustering.W..K...3.)/length(clusters_label$spectralClustering.W..K...3.)
round(SNF_Accuracy,9)

# Tumor
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0
  
# Normal
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# PVTT
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0
  
AvgF <- round((Fg1+Fg2+Fg3)/3,2)

#' 
#' 
## ------------------------------------------------------------------------
save.image(file = "HPC_SNFTool_Results_3labels.RData")

