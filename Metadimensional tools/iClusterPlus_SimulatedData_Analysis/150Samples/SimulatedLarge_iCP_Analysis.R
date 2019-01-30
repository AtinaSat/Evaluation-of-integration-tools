## iClusterPlus multi-class classification analysis for simulated large datasets group 
# We see how good the results are for 3 cluster solution. Hence k =2 as number of clusters=k+1

library(iClusterPlus)
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


setwd("/home/n9870130/Benchmark_Integrative_Tools/simulated_data/large_datasets_group/IClusterPlus_Analysis")

# load simulated data large dataset A 
load("/home/n9870130/Benchmark_Integrative_Tools/simulated_data/large_datasets_group/SimulatedCaseA_NN_150Samples.RData")
clustersLANN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LANN")


# load simulated data large dataset B
load("/home/n9870130/Benchmark_Integrative_Tools/simulated_data/large_datasets_group/SimulatedCaseB_NN_150Samples.RData")
clustersLBNN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LBNN")


# load simulated data large dataset C
load("/home/n9870130/Benchmark_Integrative_Tools/simulated_data/large_datasets_group/SimulatedCaseC_NN_150Samples.RData")
clustersLCNN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LCNN")


# load simulated data large dataset D
load("/home/n9870130/Benchmark_Integrative_Tools/simulated_data/large_datasets_group/SimulatedCaseD_NN_150Samples.RData")
clustersLDNN <- iCPAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "LDNN")


save(clustersLANN, clustersLBNN, clustersLCNN, clustersLDNN, file = "SimulatedLargeDatasetsGroup_iCP_Analysis.RData")

