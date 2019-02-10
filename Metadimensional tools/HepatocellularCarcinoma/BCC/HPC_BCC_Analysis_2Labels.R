#' ---
#' title: "Hepatocellular Carcinoma BCC - 2 labels"
#' author: "Anita"
#' date: "7 September 2018"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
## ------------------------------------------------------------------------
library(bayesCC)
library(tictoc)
load("../../HPC_filtered_data.Rdata")

#' 
## ------------------------------------------------------------------------
mrna2 <- t(log(mrna+1))
me2 <- t(me)
mirna2 <- t(log(mirna+1))
X = list(mrna2, me2, mirna2)

#' 
## ------------------------------------------------------------------------
lapply(X, dim)

#Fit model for K=2
Results = list()
set.seed(999)
tic("BCC HPC2Labels")
Results[[1]] = bayesCC(X,K=2)
toc()

K2 = Results[[1]]

Clusters = rep(1,60)
Clusters[K2$Cbest[,2]==1] = 2
Clusters <- data.frame(Clusters)
rownames(Clusters) <- data.labels$Samples

#' 
#' Performance analysis
#' 
## ------------------------------------------------------------------------
Clusters$Truth <- data.labels$Samples
Clusters$Truth <- gsub("N.{1,4}", "2", Clusters$Truth)
Clusters$Truth <- gsub("T.{1,4}", "1", Clusters$Truth)
Clusters$Truth <- gsub("P.{1,4}", "1", Clusters$Truth)

confusion.mat <- table(Truth = Clusters$Truth, Pred = Clusters$Clusters)
confusion.mat

#' 
## ------------------------------------------------------------------------
BCC_Accuracy = sum(Clusters$Truth == Clusters$Clusters)/length(Clusters$Clusters)
round(BCC_Accuracy,4)

# Tumour/PVTT
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0
  
# Normal
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

AvgF <- round((Fg1+Fg2)/2,2)
  
#save.image(file = "HPC_BCC_Results_2labels.RData")

