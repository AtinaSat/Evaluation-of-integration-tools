# BCC works only when further filtering of the data is done.
# 
## ------------------------------------------------------------------------
library(bayesCC)
load("/home/anita/Benchmarking/GBM/Data/GBM_filtered_data.Rdata")

# Data details
## ------------------------------------------------------------------------
# mrna2 is log2(affy RMA)
mrna2 <- t(gbm.dat.rna)
me2 <- t(gbm.dat.me)
protein2 <- t(gbm.dat.protein)
X = list(mrna2, me2, protein2)

# 
## ------------------------------------------------------------------------
lapply(X, dim)
#Fit model for K=4
Results = list()
set.seed(7439285)
tic("BCC GBM")
Results[[1]] = BCC(X,K=4)  
toc()
K4 = Results[[1]]

# BCC GBM: 195.293 sec elapsed

Clusters = rep(1,100)
Clusters[K4$Cbest[,2]==1] = 2
Clusters[K4$Cbest[,3]==1] = 3
Clusters[K4$Cbest[,4]==1] = 4
# 
## ------------------------------------------------------------------------
clusters_label <- data.frame(data.labels$Label)
colnames(clusters_label) <- "Samples"
clusters_label$Samples <- as.character(clusters_label$Samples)
table(clusters_label$Samples)

clusters_label$Truth <- clusters_label$Samples
clusters_label$Truth <- gsub("Classical", "3", clusters_label$Truth)
clusters_label$Truth <- gsub("Mesenchymal", "4", clusters_label$Truth)
clusters_label$Truth <- gsub("Neural", "2", clusters_label$Truth)
clusters_label$Truth <- gsub("Proneural", "1", clusters_label$Truth)

head(clusters_label)

clusters_label$Clusters <- Clusters

confusion.mat <- table(Truth = clusters_label$Truth, Pred = clusters_label$Clusters)
confusion.mat

# 
## ------------------------------------------------------------------------
BCC_Accuracy = sum(clusters_label$Truth == clusters_label$Clusters)/length(clusters_label$Clusters)
round(BCC_Accuracy,4)


# Proneural Samples (Cluster 1)
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0

# Neural (Cluster 2) NO samples in cluster 2
Pg2 <- 0 # Precision
Rg2 <- 0 # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# Classical (Cluster 3)
Pg3 <- round(confusion.mat[3,2]/sum(confusion.mat[,2]),2) # Precision
Rg3 <- round(confusion.mat[3,2]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0

# Mesenchymal samples (Cluster 4)
Pg4 <- round(confusion.mat[4,3]/sum(confusion.mat[,3]),2) # Precision
Rg4 <- round(confusion.mat[4,3]/sum(confusion.mat[4,]),2) # Recall
Fg4 <- round(2*((Pg4*Rg4)/(Pg4+Rg4)),2) # F1-score
Fg4[is.nan(Fg4)] <- 0


AvgF <- round((Fg1+Fg2+Fg3+Fg4)/4,2)
setwd("/home/anita/Benchmarking/GBM/BCC")
save.image(file = "GBM_BCC_Results.RData")

