## iClusterPlus multi-class classification analysis for simulated small datasets group 
# We see how good the results are for 3 cluster solution. Hence k =2 as number of clusters=k+1

library(iClusterPlus)
library(lattice)
library(tictoc)

# Create functions to run
iCPAnalysis <- function(geSimulated, mirnaSimulated, methSimulated, sample_labels, dataname){
  set.seed(999)
  k = 2 # number of clusters is k+1
  
  tic(dataname)
  cv.fit <- tune.iClusterPlus(dt1 = geSimulated, dt2 = mirnaSimulated,
                              dt3 = methSimulated, type = c( "gaussian","gaussian", "gaussian"),
                              n.lambda = NULL, K=k)
  toc()
  
  output = alist()
  output[[1]] = cv.fit
  
  nLambda = nrow(output[[1]]$lambda)
  nK = length(output)
  BIC = getBIC(output)
  devR = getDevR(output)

  minBICid = apply(BIC,2,which.min)
  devRatMinBIC = rep(NA,nK)
  for(i in 1:nK){
    devRatMinBIC[i] = devR[minBICid[i],i]
  }
  
  # best fit based on lambda
  best.fit=output[[1]]$fit[[which.min(BIC[,1])]]
   
  clusters_label <- data.frame(best.fit$clusters)
  colnames(clusters_label)[1] <- "Predicted"
  clusters_label$Truth <- sample_labels
  clusters_label$Truth <- gsub("Group1", "1", clusters_label$Truth)
  clusters_label$Truth <- gsub("Group2", "2", clusters_label$Truth)
  clusters_label$Truth <- gsub("Group3", "3", clusters_label$Truth)
  
  return(clusters_label)
} 


# load simulated small dataset A 
load("../../simulated_data/simulatedDatasets/SimulatedCaseA_NN.RData")
clustersANN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "ANN")


# load simulated small dataset B
load("../../simulated_data/simulatedDatasets/SimulatedCaseB_NN.RData")
clustersBNN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "BNN")


# load simulated small dataset C
load("../../simulated_data/simulatedDatasets/SimulatedCaseC_NN.RData")
clustersCNN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "CNN")


# load simulated small dataset D
load("../../simulated_data/simulatedDatasets/SimulatedCaseD_NN.RData")
clustersDNN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "DNN")


save(clustersANN, clustersBNN, clustersCNN, clustersDNN, file = "SimulatedSmallDatasetsGroup_iCP_Analysis.RData")

