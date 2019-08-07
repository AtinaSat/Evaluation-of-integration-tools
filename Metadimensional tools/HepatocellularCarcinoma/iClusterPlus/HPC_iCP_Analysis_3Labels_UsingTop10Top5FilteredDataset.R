## iClusterPlus multi-class classification analysis for HPC data Top10Top5Filtered_data
# We how good the results are for 3 cluster solution. Hence k =2 as number of clusters=k+1

library(iClusterPlus)
library(lattice)
library(tictoc)

# Loading data
load("/home/n9870130/Benchmark_Integrative_Tools/HPC_Top25Filtered/data/HPC_Top10Top5Filtered_data.Rdata")

# formatting data
mrna2 = as.matrix(mrna)
mirna2 = as.matrix(mirna)
me = as.matrix(me)

rm(mirna,mrna)

## Running iClusterPlus
set.seed(999)
k = 2 # number of clusters is k+1

tic("iClusterMultiClass")
cv.fit = tune.iClusterPlus(cpus = 24, dt1 = mrna2, dt2 = mirna2,
                             dt3 = me, type = c( "gaussian","gaussian", "gaussian"),
                             n.lambda = NULL, K=k)
toc()

output = alist()
output[[1]] <- cv.fit

nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}


best.fit=output[[1]]$fit[[which.min(BIC[,1])]]

# important features if required later
features = alist()
features[[1]] <- colnames(mrna2)
features[[2]] <- colnames(mirna2)
features[[3]] <- colnames(me)
sigfeatures = alist()
for (i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures) <- c("mRNA", "miRNA", "Methylation")

dataname = "HPC_3labels_UsingTop10Top5Filtered_data"
save(best.fit, sigfeatures, file = paste(dataname,".best.cv.fit.k", k,".RData", sep = ""))

clusters_label <- getClusters(output)
clusters_label <- data.frame(clusters_label)
colnames(clusters_label)[1] <- "Predicted"

## Rerun after loading HPC3Labels_iCPAnalysis.RData

clustersHPC3Labels = clusters_label
clustersHPC3Labels$Truth1 <- data.labels$Samples
clustersHPC3Labels$Truth1 <- gsub("N.{1,4}", "N", clustersHPC3Labels$Truth1)
clustersHPC3Labels$Truth1 <- gsub("T.{1,4}", "T", clustersHPC3Labels$Truth1)
clustersHPC3Labels$Truth1 <- gsub("P.{1,4}", "P", clustersHPC3Labels$Truth1)
confusion.mat =table(truth=clustersHPC3Labels$Truth1, predicted = clustersHPC3Labels$Predicted)
confusion.mat

#       predicted
# truth  1  2  3
#     N  1  0 19
#     P  8  9  3
#     T  9  9  2

## Cluster 1: Tumour Samples; Cluster 2: PVTT Samples; Cluster 3: Normal Samples


clustersHPC3Labels = clusters_label
clustersHPC3Labels$Truth <- data.labels$Samples
clustersHPC3Labels$Truth <- gsub("N.{1,4}", "3", clustersHPC3Labels$Truth)
clustersHPC3Labels$Truth <- gsub("T.{1,4}", "1", clustersHPC3Labels$Truth)
clustersHPC3Labels$Truth <- gsub("P.{1,4}", "2", clustersHPC3Labels$Truth)
confusion.mat =table(truth=clustersHPC3Labels$Truth, predicted = clustersHPC3Labels$Predicted)
confusion.mat

#       predicted
# truth  1  2  3
#     1  9  9  2
#     2  8  9  3
#     3  1  0 19


# Tumour samples (Cluster 1)
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0
  

# PVTT samples (Cluster 2)
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

# Normal Samples (Cluster 3)
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0
  
Accuracy <-  sum(diag(confusion.mat))/60
Avg.F1 <- (Fg1+Fg2+Fg3)/3

save.image(file = "HPC3Labels_iCPAnalysis_UsingTop10Top5Filtered_data.RData")
