#' ---
#' title: "SNF Analysis using Simulated Data (150 Samples)"
#' author: "Anita"
#' date: "31 August 2018"
#' output:
#'   pdf_document: default
#'   html_document: default
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # Introduction
#' Simulated dataset was used to anaylse the classifying power of Similarity Network Fusion (SNF). The dataset was simulated from breast cancer data from TCGA.  
#'   
#' In this dataset there are three data types, gene expression, miRNA expression, and methylation, each simulated as normal distribution from real data. There are 60 samples and three data labels, "Group1", "Group2", "Group3".   
#'   
#' Gene expression, miRNA expression and methylation have 500 features, each. In case of noisy data, 100 noisy features have been added to gene expression, and 20 noisy features to miRNA and methylation data, respectively. Each case has raw data and raw data + noisy data.  
#' 
#' ##### Simulated dataset case A 
#' All three groups are clearly distinguishable  
#' 
#' ##### Simulated dataset case B
#' Data type 1 (Gene expression): Group3 has mixed signal  
#' Data type 2 (miRNA): Group1 has mixed signal  
#' Data type 3 (methylation): All three groups are clearly distinguishable  
#' 
#' ##### Simulated dataset case C
#' Data type 1 (Gene expression): Group3 has mixed signal  
#' Data type 2 (miRNA): Group2 has mixed signal  
#' Data type 3 (methylation): Group2 has mixed signal  
#' 
#' ##### Simulated dataset case D
#' Data type 1 (Gene expression): Group3 has mixed signal    
#' Data type 2 (miRNA): Group2 has mixed signal  
#' Data type 3 (methylation): Group1 has mixed signal  
#' 
## ------------------------------------------------------------------------
rm(list = ls())
library(SNFtool)
library(tictoc)

#' 
#' 
## ------------------------------------------------------------------------
SNFAnalysis <- function(geSimulated, mirnaSimulated, methSimulated, sample_labels, dataname){
  
  # calculate distance
  dist1 <- as.matrix(dist(geSimulated))
  dist2 <- as.matrix(dist(mirnaSimulated))
  dist3 <- as.matrix(dist(methSimulated))
  
  tic(dataname)
  # calculate affinity using default values
  W1 <- affinityMatrix(dist1)
  W2 <- affinityMatrix(dist2)
  W3 <- affinityMatrix(dist3)
  
  # combine clusters from each data type using default values
  W = SNF(list(W1,W2, W3))
  
  # spectral clustering default type
  clusters_label <- data.frame(spectralClustering(W,K = 3))
  toc()
  
  colnames(clusters_label)[1] <- "Predicted"
  clusters_label$Truth <- sample_labels
  clusters_label$Truth <- gsub("Group1", "1", clusters_label$Truth)
  clusters_label$Truth <- gsub("Group2", "2", clusters_label$Truth)
  clusters_label$Truth <- gsub("Group3", "3", clusters_label$Truth)
  
  return(clusters_label)
} 

PerformanceAnalysis <- function(confusion.mat){
  # Calculating the precision, recall, and F scores
  # Group1
  Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
  Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
  Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
  
  g1 <- c(Pg1, Rg1, Fg1)

  # Group2
  Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
  Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
  Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score

  g2 <- c(Pg2, Rg2, Fg2)

  # Group3
  Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
  Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
  Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
  
  g3 <- c(Pg3, Rg3, Fg3)

  AvgF <- round((Fg1+Fg2+Fg3)/3,2)
  
  res <- list(g1, g2, g3, AvgF)
  
  return(res)
}


#' 
#' ### Large datasets group
#' #### Simulated dataset A
#' 
## ------------------------------------------------------------------------
# load simulated dataset A no noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseA_NN_150Samples.RData")

clustersANN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LANN")
confusion.mat = table(Truth = clustersANN$Truth, Pred = clustersANN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 2, Cluster 2 = Group 3, Cluster 3 = Group 1
# Renaming groups with cluster numbers
predlabels = clustersANN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group2", "Group3", "Group1")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersANN$Truth <- sample_labels
confusion.mat = table(Truth = clustersANN$Truth, predlabels)
confusion.mat

AccuracyANN = round(sum(clustersANN$Truth == predlabels)/length(predlabels),2)
resANN <- PerformanceAnalysis(confusion.mat)

GroupScoresSNF <- data.frame(matrix(,nrow = 24, ncol = 6))
colnames(GroupScoresSNF) <- c("Tool", "Dataset", "Group", "Precision", "Recall", "FScore")
GroupScoresSNF[,1] <- "SNF"
GroupScoresSNF[,3] <- rep(c("Group1", "Group2", "Group3"), 8)

DataScoresSNF <- data.frame(matrix(,nrow = 8, ncol = 4))
colnames(DataScoresSNF) <- c("Tool", "Dataset", "Accuracy", "AvgFScore")
DataScoresSNF[,1] <- "SNF"

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[1:3,2] <- "A_NN"
GroupScoresSNF[1,4:6] <- resANN[[1]]
GroupScoresSNF[2,4:6] <- resANN[[2]]
GroupScoresSNF[3,4:6] <- resANN[[3]]

DataScoresSNF[1,2] <- "A_NN"
DataScoresSNF[1,3] <- AccuracyANN
DataScoresSNF[1,4] <- resANN[[4]]

#' 
#' ### Large noisy datasets group
#' #### Simulated dataset A
#' 
## ------------------------------------------------------------------------
# load simulated dataset A with noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseA_WN_150Samples.RData")

clustersAWN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LAWN")
confusion.mat = table(Truth = clustersAWN$Truth, Pred = clustersAWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 2, Cluster 2 = Group 3, Cluster 3 = Group 1
# Renaming groups with cluster numbers
predlabels = clustersAWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group2", "Group3", "Group1")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersAWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersAWN$Truth, predlabels)
confusion.mat

AccuracyAWN = round(sum(clustersAWN$Truth == predlabels)/length(predlabels),2)
resAWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[4:6,2] <- "A_WN"
GroupScoresSNF[4,4:6] <- resAWN[[1]]
GroupScoresSNF[5,4:6] <- resAWN[[2]]
GroupScoresSNF[6,4:6] <- resAWN[[3]]

DataScoresSNF[2,2] <- "A_WN"
DataScoresSNF[2,3] <- AccuracyAWN
DataScoresSNF[2,4] <- resAWN[[4]]

#' 
#' ### Large datasets group
#' #### Simulated dataset B
#' 
## ------------------------------------------------------------------------
# load simulated data case B no noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseB_NN_150Samples.RData")

clustersBNN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LBNN")
confusion.mat = table(Truth = clustersBNN$Truth, Pred = clustersBNN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 1, Cluster 2 = Group 3, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersBNN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group1", "Group3", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersBNN$Truth <- sample_labels
confusion.mat = table(Truth = clustersBNN$Truth, predlabels)
confusion.mat

AccuracyBNN = round(sum(clustersBNN$Truth == predlabels)/length(predlabels),2)
resBNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[7:9,2] <- "B_NN"
GroupScoresSNF[7,4:6] <- resBNN[[1]]
GroupScoresSNF[8,4:6] <- resBNN[[2]]
GroupScoresSNF[9,4:6] <- resBNN[[3]]

DataScoresSNF[3,2] <- "B_NN"
DataScoresSNF[3,3] <- AccuracyBNN
DataScoresSNF[3,4] <- resBNN[[4]]

#' 
#' ### Large noisy datasets group
#' #### Simulated dataset B
#' 
## ------------------------------------------------------------------------
# load simulated data case B with noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseB_WN_150Samples.RData")

clustersBWN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LBWN")
confusion.mat = table(Truth = clustersBWN$Truth, Pred = clustersBWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 1, Cluster 2 = Group 3, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersBWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group1", "Group3", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersBWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersBWN$Truth, predlabels)
confusion.mat

AccuracyBWN = round(sum(clustersBWN$Truth == predlabels)/length(predlabels),2)
resBWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[10:12,2] <- "B_WN"
GroupScoresSNF[10,4:6] <- resBWN[[1]]
GroupScoresSNF[11,4:6] <- resBWN[[2]]
GroupScoresSNF[12,4:6] <- resBWN[[3]]

DataScoresSNF[4,2] <- "B_WN"
DataScoresSNF[4,3] <- AccuracyBWN
DataScoresSNF[4,4] <- resBWN[[4]]

#' 
#' ### Large datasets group
#' #### Simulated dataset C
#' 
## ------------------------------------------------------------------------
# load simulated data case C no noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseC_NN_150Samples.RData")

clustersCNN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LCNN")
confusion.mat = table(Truth = clustersCNN$Truth, Pred = clustersCNN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 1, Cluster 2 = Group 2, Cluster 3 = Group 3
# Renaming groups with cluster numbers
predlabels = clustersCNN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group1", "Group2", "Group3")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersCNN$Truth <- sample_labels
confusion.mat = table(Truth = clustersCNN$Truth, predlabels)
confusion.mat

AccuracyCNN = round(sum(clustersCNN$Truth == predlabels)/length(predlabels),2)
resCNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[13:15,2] <- "C_NN"
GroupScoresSNF[13,4:6] <- resCNN[[1]]
GroupScoresSNF[14,4:6] <- resCNN[[2]]
GroupScoresSNF[15,4:6] <- resCNN[[3]]

DataScoresSNF[5,2] <- "C_NN"
DataScoresSNF[5,3] <- AccuracyCNN
DataScoresSNF[5,4] <- resCNN[[4]]

#' 
#' ### Large noisy datasets group
#' #### Simulated dataset C
#' 
## ------------------------------------------------------------------------
# load simulated data case C with noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseC_WN_150Samples.RData")

clustersCWN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LCWN")
confusion.mat = table(Truth = clustersCWN$Truth, Pred = clustersCWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 2, Cluster 2 = Group 3, Cluster 3 = Group 1
# Renaming groups with cluster numbers
predlabels = clustersCWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group2", "Group3", "Group1")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersCWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersCWN$Truth, predlabels)
confusion.mat

AccuracyCWN = round(sum(clustersCWN$Truth == predlabels)/length(predlabels),2)
resCWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[16:18,2] <- "C_WN"
GroupScoresSNF[16,4:6] <- resCWN[[1]]
GroupScoresSNF[17,4:6] <- resCWN[[2]]
GroupScoresSNF[18,4:6] <- resCWN[[3]]

DataScoresSNF[6,2] <- "C_WN"
DataScoresSNF[6,3] <- AccuracyCWN
DataScoresSNF[6,4] <- resCWN[[4]]

#' 
#' ### Large datasets group
#' #### Simulated dataset D
#' 
## ------------------------------------------------------------------------
# load simulated data case D no noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseD_NN_150Samples.RData")

clustersDNN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LDNN")
confusion.mat = table(Truth = clustersDNN$Truth, Pred = clustersDNN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 1, Cluster 2 = Group 3, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersDNN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group1", "Group3", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersDNN$Truth <- sample_labels
confusion.mat = table(Truth = clustersDNN$Truth, predlabels)
confusion.mat

AccuracyDNN = round(sum(clustersDNN$Truth == predlabels)/length(predlabels),2)
resDNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[19:21,2] <- "D_NN"
GroupScoresSNF[19,4:6] <- resDNN[[1]]
GroupScoresSNF[20,4:6] <- resDNN[[2]]
GroupScoresSNF[21,4:6] <- resDNN[[3]]

DataScoresSNF[7,2] <- "D_NN"
DataScoresSNF[7,3] <- AccuracyDNN
DataScoresSNF[7,4] <- resDNN[[4]]

#' 
#' ### Large noisy datasets group
#' #### Simulated dataset D
#' 
## ------------------------------------------------------------------------
# load simulated data case D with noise
load("/home/anita/Benchmarking/simulated_data/SimulatedCaseD_WN_150Samples.RData")

clustersDWN <- SNFAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LDWN")
confusion.mat = table(Truth = clustersDWN$Truth, Pred = clustersDWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 1, Cluster 2 = Group 2, Cluster 3 = Group 3
# Renaming groups with cluster numbers
predlabels = clustersDWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group1", "Group2", "Group3")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersDWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersDWN$Truth, predlabels)
confusion.mat

AccuracyDWN = round(sum(clustersDWN$Truth == predlabels)/length(predlabels),2)
resDWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresSNF[22:24,2] <- "D_WN"
GroupScoresSNF[22,4:6] <- resDWN[[1]]
GroupScoresSNF[23,4:6] <- resDWN[[2]]
GroupScoresSNF[24,4:6] <- resDWN[[3]]

DataScoresSNF[8,2] <- "D_WN"
DataScoresSNF[8,3] <- AccuracyDWN
DataScoresSNF[8,4] <- resDWN[[4]]

#' 
#' ### Plots
#' 
## ------------------------------------------------------------------------
library(ggplot2)

DataScoresSNF$Data <- rep(c("A", "B", "C", "D"), each = 2)
DataScoresSNF$Noise <- rep(c("NN", "WN"), 4)

p <- ggplot(DataScoresSNF, aes(x = Data, y = AvgFScore, fill = Noise)) + 
  geom_bar(stat = "identity", position = position_dodge())+ 
  labs(title="SNF",x="Dataset", y = "Average F1 Score")+
  ylim(0,1)+
  theme(text = element_text(size=20), axis.text = element_text(face = "bold", color = "black"))+
  geom_text(aes(label=AvgFScore), position=position_dodge(width=0.9), vjust=-0.25)
p

p2 <- ggplot(DataScoresSNF, aes(x = Data, y = Accuracy, fill = Noise)) + 
  geom_bar(stat = "identity", position = position_dodge())+ 
  labs(title="SNF",x="Dataset", y = "Accuracy")+
  ylim(0,1)+
  theme(text = element_text(size=20), axis.text = element_text(face = "bold", color = "black"))+
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.25)
p2

GroupScoresSNF$Data <- rep(c("A", "B", "C", "D"), each = 6)
GroupScoresSNF$Noise <- rep(rep(c("NN", "WN"), each = 3),4)
GroupScoresSNF$Group <- factor(GroupScoresSNF$Group, levels = c('Group1', 'Group2', 'Group3'))

p3 <-ggplot(GroupScoresSNF, aes(x = Dataset, y = FScore, fill = Group)) + 
  geom_bar(stat = "identity", position = position_dodge())+ 
  labs(title="SNF",x="", y = "F1 Score")+
  ylim(0,1)+
  theme(text = element_text(size=14), axis.text = element_text(face = "bold", color = "black"))+
  geom_text(aes(label=FScore), position=position_dodge(width=0.9), vjust=-0.25, size = 4)
p3 

#' 
## ------------------------------------------------------------------------
setwd("/home/anita/Benchmarking/simulated_data/SNF_Simulated_Analysis/150Samples/")
write.table(DataScoresSNF, file = "DataScoresSNF_150Samples.txt", quote = F, sep = "\t")
write.table(GroupScoresSNF, file = "GroupScoresSNF_150Samples.txt", quote = F, sep = "\t")
save.image(file = "SNF_Simulated_Analysis_Results_150Samples.RData")

