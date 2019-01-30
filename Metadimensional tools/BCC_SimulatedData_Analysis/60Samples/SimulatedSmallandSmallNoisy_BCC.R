#' ---
#' title: "BCC Analysis using Simulated Data"
#' author: "Anita"
#' date: "20 August 2018"
#' output: html_document
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' 
#' # Introduction
#' Simulated dataset was used to anaylse the classifying power of BCC. The dataset was simulated from breast cancer data from TCGA.  
#'   
#' In this dataset there are three data types, gene expression, miRNA expression, and methylation, each simulated as normal distribution from real data. There are 60 samples and three data labels, "Group1", "Group2", "Group3".   
#'   
#' Gene expression, miRNA expression and methylation have 500 features, each. In case of noisy data, 100 noisy features have been added to gene expression, and 20 noisy features to miRNA and methylation data, respectively. Each case has raw data and raw data + noisy data.  
#' 
#' ##### Simulated dataset A 
#' All three groups are clearly distinguishable  
#' 
#' ##### Simulated dataset B
#' Data type 1 (Gene expression): Group3 has mixed signal  
#' Data type 2 (miRNA): Group1 has mixed signal  
#' Data type 3 (methylation): All three groups are clearly distinguishable  
#' 
#' ##### Simulated dataset C
#' Data type 1 (Gene expression): Group3 has mixed signal  
#' Data type 2 (miRNA): Group2 has mixed signal  
#' Data type 3 (methylation): Group2 has mixed signal  
#' 
#' ##### Simulated dataset D
#' Data type 1 (Gene expression): Group3 has mixed signal    
#' Data type 2 (miRNA): Group2 has mixed signal  
#' Data type 3 (methylation): Group1 has mixed signal  
#' 
## ------------------------------------------------------------------------
library(gtools)
library(gplots)
library(bayesCC)
rm(list = ls())

#' 
## ------------------------------------------------------------------------
BCCAnalysis <- function(geSimulated, mirnaSimulated, methSimulated, sample_labels, dataname){
  
  ge2 <- t(geSimulated)
  mirna2 <- t(mirnaSimulated)
  me2 <- t(methSimulated)
  
  X = list(ge2, mirna2, me2)
  
  #Fit model for K=3 
  Results = list()
  set.seed(123)
  
  tic(dataname)
  Results[[1]] = bayesCC(X,K=3)
  toc()
  
  K3 = Results[[1]]
  
  clusters_label = rep(1,60)
  clusters_label[K3$Cbest[,2]==1] = 2
  clusters_label[K3$Cbest[,3]==1] = 3
  
  clusters_label <- data.frame(clusters_label)
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
  Fg1[is.nan(Fg1)] <- 0
  
  g1 <- c(Pg1, Rg1, Fg1)

  # Group2
  Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
  Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
  Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
  Fg2[is.nan(Fg2)] <- 0

  g2 <- c(Pg2, Rg2, Fg2)

  # Group3
  Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
  Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
  Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
  Fg3[is.nan(Fg3)] <- 0
  
  g3 <- c(Pg3, Rg3, Fg3)

  AvgF <- round((Fg1+Fg2+Fg3)/3,2)
  
  res <- list(g1, g2, g3, AvgF)
  
  return(res)
}


#' 
#' ### Small datasets group
#' #### Simulated dataset A 
#' 
## ------------------------------------------------------------------------
# load simulated data case A no noise
Evaluation-of-integration-tools/Metadimensional tools/simulated_data/simulatedDatasets/
load("../../simulated_data/simulatedDatasets/SimulatedCaseA_NN.RData")

clustersANN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "ANN")

confusion.mat = table(Truth = clustersANN$Truth, Pred = clustersANN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 2, Cluster 2 = Group 1, Cluster 3 = Group 3
# Renaming groups with cluster numbers
predlabels = clustersANN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group2", "Group1", "Group3")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersANN$Truth <- sample_labels
confusion.mat = table(Truth = clustersANN$Truth, predlabels)
confusion.mat

AccuracyANN = round(sum(clustersANN$Truth == predlabels)/length(predlabels),2)
resANN <- PerformanceAnalysis(confusion.mat)

GroupScoresBCC <- data.frame(matrix(,nrow = 24, ncol = 6))
colnames(GroupScoresBCC) <- c("Tool", "Dataset", "Group", "Precision", "Recall", "FScore")
GroupScoresBCC[,1] <- "BCC"
GroupScoresBCC[,3] <- rep(c("Group1", "Group2", "Group3"), 8)

DataScoresBCC <- data.frame(matrix(,nrow = 8, ncol = 4))
colnames(DataScoresBCC) <- c("Tool", "Dataset", "Accuracy", "AvgFScore")
DataScoresBCC[,1] <- "BCC"

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[1:3,2] <- "A_NN"
GroupScoresBCC[1,4:6] <- resANN[[1]]
GroupScoresBCC[2,4:6] <- resANN[[2]]
GroupScoresBCC[3,4:6] <- resANN[[3]]

DataScoresBCC[1,2] <- "A_NN"
DataScoresBCC[1,3] <- AccuracyANN
DataScoresBCC[1,4] <- resANN[[4]]

#' 
#' ### Small noisy datasets group
#' #### Simulated dataset A 
#' 
## ------------------------------------------------------------------------
# load simulated data A with noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseA_WN.RData")


clustersAWN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "AWN")

confusion.mat = table(Truth = clustersAWN$Truth, Pred = clustersAWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 3, Cluster 2 = Group 1, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersAWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group3", "Group1", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersAWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersAWN$Truth, predlabels)
confusion.mat

AccuracyAWN = round(sum(clustersAWN$Truth == predlabels)/length(predlabels),2)
resAWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[4:6,2] <- "A_WN"
GroupScoresBCC[4,4:6] <- resAWN[[1]]
GroupScoresBCC[5,4:6] <- resAWN[[2]]
GroupScoresBCC[6,4:6] <- resAWN[[3]]

DataScoresBCC[2,2] <- "A_WN"
DataScoresBCC[2,3] <- AccuracyAWN
DataScoresBCC[2,4] <- resAWN[[4]]

#' 
#' ### Small datasets group
#' #### Simulated dataset B 
#' 
## ------------------------------------------------------------------------
# load simulated dataset B no noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseB_NN.RData")


clustersBNN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "BNN")

confusion.mat = table(Truth = clustersBNN$Truth, Pred = clustersBNN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 2, Cluster 2 = Group 1, Cluster 3 = Group 3
# Renaming groups with cluster numbers
predlabels = clustersBNN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group2", "Group1", "Group3")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersBNN$Truth <- sample_labels
confusion.mat = table(Truth = clustersBNN$Truth, predlabels)
confusion.mat

AccuracyBNN = round(sum(clustersBNN$Truth == predlabels)/length(predlabels),2)
resBNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[7:9,2] <- "B_NN"
GroupScoresBCC[7,4:6] <- resBNN[[1]]
GroupScoresBCC[8,4:6] <- resBNN[[2]]
GroupScoresBCC[9,4:6] <- resBNN[[3]]

DataScoresBCC[3,2] <- "B_NN"
DataScoresBCC[3,3] <- AccuracyBNN
DataScoresBCC[3,4] <- resBNN[[4]]

#' 
#' ### Small noisy datasets group
#' #### Simulated dataset B 
#' 
## ------------------------------------------------------------------------
# load simulated dataset B with noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseB_WN.RData")


clustersBWN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "BWN")

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
GroupScoresBCC[10:12,2] <- "B_WN"
GroupScoresBCC[10,4:6] <- resBWN[[1]]
GroupScoresBCC[11,4:6] <- resBWN[[2]]
GroupScoresBCC[12,4:6] <- resBWN[[3]]

DataScoresBCC[4,2] <- "B_WN"
DataScoresBCC[4,3] <- AccuracyBWN
DataScoresBCC[4,4] <- resBWN[[4]]

#' 
#' ### Small datasets group
#' #### Simulated dataset C
#' 
## ------------------------------------------------------------------------
# load simulated dataset C no noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseC_NN.RData")


clustersCNN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "CNN")

confusion.mat = table(Truth = clustersCNN$Truth, Pred = clustersCNN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 2, Cluster 2 = Group 3, Cluster 3 = Group 1
# Renaming groups with cluster numbers
predlabels = clustersCNN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group2", "Group1", "Group3")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersCNN$Truth <- sample_labels
confusion.mat = table(Truth = clustersCNN$Truth, predlabels)
confusion.mat

AccuracyCNN = round(sum(clustersCNN$Truth == predlabels)/length(predlabels),2)
resCNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[13:15,2] <- "C_NN"
GroupScoresBCC[13,4:6] <- resCNN[[1]]
GroupScoresBCC[14,4:6] <- resCNN[[2]]
GroupScoresBCC[15,4:6] <- resCNN[[3]]

DataScoresBCC[5,2] <- "C_NN"
DataScoresBCC[5,3] <- AccuracyCNN
DataScoresBCC[5,4] <- resCNN[[4]]

#' 
#' ### Small noisy datasets group
#' #### Simulated dataset C 
#' 
## ------------------------------------------------------------------------
# load simulated dataset C with noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseC_WN.RData")

clustersCWN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "CWN")

confusion.mat = table(Truth = clustersCWN$Truth, Pred = clustersCWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 1, Cluster 2 = Group 3, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersCWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group1", "Group3", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersCWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersCWN$Truth, predlabels)
confusion.mat

AccuracyCWN = round(sum(clustersCWN$Truth == predlabels)/length(predlabels),2)
resCWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[16:18,2] <- "C_WN"
GroupScoresBCC[16,4:6] <- resCWN[[1]]
GroupScoresBCC[17,4:6] <- resCWN[[2]]
GroupScoresBCC[18,4:6] <- resCWN[[3]]

DataScoresBCC[6,2] <- "C_WN"
DataScoresBCC[6,3] <- AccuracyCWN
DataScoresBCC[6,4] <- resCWN[[4]]

#' 
#' ### Small datasets group
#' #### Simulated dataset D
#' 
## ------------------------------------------------------------------------
# load simulated dataset D no noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseD_NN.RData")


clustersDNN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "DNN")

confusion.mat = table(Truth = clustersDNN$Truth, Pred = clustersDNN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 3, Cluster 2 = Group 1, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersDNN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group3", "Group1", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))


clustersDNN$Truth <- sample_labels

confusion.mat = table(Truth = clustersDNN$Truth, predlabels)
confusion.mat

AccuracyDNN = round(sum(clustersDNN$Truth == predlabels)/length(predlabels),2)
resDNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[19:21,2] <- "D_NN"
GroupScoresBCC[19,4:6] <- resDNN[[1]]
GroupScoresBCC[20,4:6] <- resDNN[[2]]
GroupScoresBCC[21,4:6] <- resDNN[[3]]

DataScoresBCC[7,2] <- "D_NN"
DataScoresBCC[7,3] <- AccuracyDNN
DataScoresBCC[7,4] <- resDNN[[4]]

#' 
#' ### Small noisy datasets group
#' #### Simulated dataset D
#' 
## ------------------------------------------------------------------------
# load simulated dataset D with noise
load("../../simulated_data/simulatedDatasets/SimulatedCaseD_WN.RData")


clustersDWN <- BCCAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "DWN")

confusion.mat = table(Truth = clustersDWN$Truth, Pred = clustersDWN$Predicted)
confusion.mat

#' 
## ------------------------------------------------------------------------
# Majority of cluster 1 = Group 3, Cluster 2 = Group 1, Cluster 3 = Group 2
# Renaming groups with cluster numbers
predlabels = clustersDWN$Predicted
predlabels <- as.factor(predlabels)
levels(predlabels) <- c("Group3", "Group1", "Group2")
predlabels <- factor(predlabels, c("Group1", "Group2", "Group3"))

clustersDWN$Truth <- sample_labels
confusion.mat = table(Truth = clustersDWN$Truth, predlabels)
confusion.mat

AccuracyDWN = round(sum(clustersDWN$Truth == predlabels)/length(predlabels),2)
resDWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresBCC[22:24,2] <- "D_WN"
GroupScoresBCC[22,4:6] <- resDWN[[1]]
GroupScoresBCC[23,4:6] <- resDWN[[2]]
GroupScoresBCC[24,4:6] <- resDWN[[3]]

DataScoresBCC[8,2] <- "D_WN"
DataScoresBCC[8,3] <- AccuracyDWN
DataScoresBCC[8,4] <- resDWN[[4]]

#' 
#' ### Plots
#' 
## ------------------------------------------------------------------------
library(ggplot2)

DataScoresBCC$Data <- rep(c("A", "B", "C", "D"), each = 2)
DataScoresBCC$Noise <- rep(c("NN", "WN"), 4)

p <- ggplot(DataScoresBCC, aes(x = Data, y = AvgFScore, fill = Noise)) + 
  geom_bar(stat = "identity", position = position_dodge())+ 
  labs(title="BCC",x="Dataset", y = "Average F1 Score")+
  ylim(0,1)+
  theme(text = element_text(size=20), axis.text = element_text(face = "bold", color = "black"))+
  geom_text(aes(label=AvgFScore), position=position_dodge(width=0.9), vjust=-0.25)
p

p2 <- ggplot(DataScoresBCC, aes(x = Data, y = Accuracy, fill = Noise)) + 
  geom_bar(stat = "identity", position = position_dodge())+ 
  labs(title="BCC",x="Dataset", y = "Accuracy")+
  ylim(0,1)+
  theme(text = element_text(size=20), axis.text = element_text(face = "bold", color = "black"))+
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.25)
p2

GroupScoresBCC$Data <- rep(c("A", "B", "C", "D"), each = 6)
GroupScoresBCC$Noise <- rep(rep(c("NN", "WN"), each = 3),4)
GroupScoresBCC$Group <- factor(GroupScoresBCC$Group, levels = c('Group1', 'Group2', 'Group3'))

p3 <-ggplot(GroupScoresBCC, aes(x = Dataset, y = FScore, fill = Group)) + 
  geom_bar(stat = "identity", position = position_dodge())+ 
  labs(title="BCC",x="", y = "F1 Score")+
  ylim(0,1)+
  theme(text = element_text(size=14), axis.text = element_text(face = "bold", color = "black"))+
  geom_text(aes(label=FScore), position=position_dodge(width=0.9), vjust=-0.25, size = 4)
p3 

#' 
#' 
## ------------------------------------------------------------------------
write.table(DataScoresBCC, file = "DataScoresBCC.txt", quote = F, sep = "\t")
write.table(GroupScoresBCC, file = "GroupScoresBCC.txt", quote = F, sep = "\t")
save.image(file = "BCC_Simulated_Analysis_Results.RData")

