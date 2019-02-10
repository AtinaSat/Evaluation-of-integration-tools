## iClusterPlus binary-class classification analysis for HPC data 
# We how good the results are for 2 cluster solution. Hence k =1 as number of clusters=k+1

library(iClusterPlus)
library(lattice)
library(tictoc)

# Loading data
load("../HPC_filtered_data.Rdata")

# formatting data
mrna = as.matrix(mrna)
mrna2 = log(mrna+1)

mirna = as.matrix(mirna)
mirna2 = log(mirna+1)

me = as.matrix(me)

rm(mirna,mrna)

## Running iClusterPlus
set.seed(999)
k = 1 # number of clusters is k+1

tic("iClusterBinaryClass")
cv.fit = tune.iClusterPlus(cpus = 16, dt1 = mrna2, dt2 = mirna2,
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

#dataname = "HPC_2labels"
#save(best.fit, sigfeatures, file = paste(dataname,".best.cv.fit.k", k,".RData", sep = ""))

clusters_label <- getClusters(output)
clusters_label <- data.frame(clusters_label)
colnames(clusters_label)[1] <- "Predicted"

clustersHPC2Labels = clusters_label
clustersHPC2Labels$Truth <- data.labels$Samples

clustersHPC2Labels$Truth <- gsub("N.{1,4}", "1", clustersHPC2Labels$Truth)
clustersHPC2Labels$Truth <- gsub("T.{1,4}", "2", clustersHPC2Labels$Truth)
clustersHPC2Labels$Truth <- gsub("P.{1,4}", "2", clustersHPC2Labels$Truth)
confusion.mat =table(truth=clustersHPC2Labels$Truth, predicted = clustersHPC2Labels$Predicted)
confusion.mat

# predicted
# truth  1  2
#     1 20  0
#     2  1 39

# Normal
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0
  
  
# PVTT/Tumour
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

AvgF <- round((Fg1+Fg2)/2,2)

Accuracy <-  sum(diag(confusion.mat))/60

save.image(file = "HPC2Labels_iCPAnalysis.RData")
