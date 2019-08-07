## HPC 2 labels SNFTool analysis using Top10Top5filtered

library(SNFtool)
library(tictoc)
load("HPC_Top10Top5Filtered_data.Rdata")

tic("SNF-HPC2Labels")
# Calculating distance matrices
dist1 <- as.matrix(dist(mrna))
dist2 <- as.matrix(dist(mirna))
dist3 <- as.matrix(dist(me))

# Calculating affinity matrics
W1 <- affinityMatrix(dist1)
W2 <- affinityMatrix(dist2)
W3 <- affinityMatrix(dist3)

# Combining the clusters
W = SNF(list(W1,W2, W3))

# Spectral clustering
clusters_label <- data.frame(spectralClustering(W,K = 2))

toc()
# 

clusters_label$Samples <- data.labels$Samples
clusters_label$Samples <- as.character(clusters_label$Samples)

clusters_label$Truth <- clusters_label$Samples
clusters_label$Truth1 <- clusters_label$Samples

clusters_label$Truth1 <- gsub("N.{1,4}", "N", clusters_label$Truth1)
clusters_label$Truth1 <- gsub("T.{1,4}", "TP", clusters_label$Truth1)
clusters_label$Truth1 <- gsub("P.{1,4}", "TP", clusters_label$Truth1)

confusion.mat = table(true = clusters_label$Truth1, Pred = clusters_label$spectralClustering.W..K...2.)
confusion.mat

#     Pred
# true  1  2
# N     0 20
# TP   37  3


# Cluster 1: Tumour/PVTT; Cluster 2: Normal; 
clusters_label$Truth <- clusters_label$Samples
clusters_label$Truth <- gsub("N.{1,4}", "2", clusters_label$Truth)
clusters_label$Truth <- gsub("T.{1,4}", "1", clusters_label$Truth)
clusters_label$Truth <- gsub("P.{1,4}", "1", clusters_label$Truth)

confusion.mat = table(true = clusters_label$Truth, Pred = clusters_label$spectralClustering.W..K...2.)
confusion.mat

#     Pred
# true  1  2
# 1    37  3
# 2     0 20

# Performance analysis
SNF_Accuracy = sum(clusters_label$Truth == clusters_label$spectralClustering.W..K...2.)/length(clusters_label$spectralClustering.W..K...2.)
round(SNF_Accuracy,9)
# 0.95

# Tumor/PVTT
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

save.image(file = "HPC_SNFTool_Results_2labels_UsingTop10Top5FilteredDataset.RData")

