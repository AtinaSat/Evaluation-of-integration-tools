# Hepatocellular carcinoma BCC analysis 3 labels using Top10Top5Filtered_HPCData

library(bayesCC)
load("HPC_Top10Top5Filtered_data.Rdata")
mrna2 <- t(mrna)
me2 <- t(me)
mirna2 <- t(mirna)

X = list(mrna2, me2, mirna2)

lapply(X, dim)
# [[1]]
# [1] 1728   60
# 
# [[2]]
# [1] 24272    60
# 
# [[3]]
# [1] 60 60

#Fit model for K=3
Results = list()
set.seed(3566)
tic("BCC HPC3Labels")
Results[[1]] = bayesCC(X,K=3)
toc()
# BCC HPC3Labels: 2396.732 sec elapsed
K3 = Results[[1]]

Clusters = rep(1,60)
Clusters[K3$Cbest[,2]==1] = 2
Clusters[K3$Cbest[,3]==1] = 3
Clusters <- data.frame(Clusters)
rownames(Clusters) <- data.labels$Samples

Clusters$Truth1 <- data.labels$Samples
Clusters$Truth1 <- gsub("N.{1,4}", "N", Clusters$Truth1)
Clusters$Truth1 <- gsub("T.{1,4}", "T", Clusters$Truth1)
Clusters$Truth1 <- gsub("P.{1,4}", "P", Clusters$Truth1)
confusion.mat <- table(Truth = Clusters$Truth1, Pred = Clusters$Clusters)
confusion.mat

#       Pred
# Truth  1  2  3
#   N    0 20  0
#   P    6  1 13
#   T    6  1 13

# Cluster 1: Either PVTT/Tumour; Cluster 2: Normal Samples; Cluster 3: PVTT/Tumour 

Clusters$Truth <- data.labels$Samples
Clusters$Truth <- gsub("N.{1,4}", "2", Clusters$Truth)
Clusters$Truth <- gsub("T.{1,4}", "3", Clusters$Truth)
Clusters$Truth <- gsub("P.{1,4}", "1", Clusters$Truth)

confusion.mat <- table(Truth = Clusters$Truth, Pred = Clusters$Clusters)
confusion.mat

#       Pred
# Truth  1  2  3
#   1    6  1 13
#   2    0 20  0
#   3    6  1 13

BCC_Accuracy = sum(Clusters$Truth == Clusters$Clusters)/length(Clusters$Clusters)
round(BCC_Accuracy,4)
# [1] 0.65

# PVTT Samples (Cluster 1)
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0
  
# Normal (Cluster 2)
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# Tumour samples (Cluster 3)
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0
  
AvgF <- round((Fg1+Fg2+Fg3)/3,2)

save.image(file = "HPC_BCC_Results_3labels_UsingTop10Top5FilteredDataset.RData")