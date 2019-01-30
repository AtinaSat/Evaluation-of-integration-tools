rm(list = ls())
library(mixOmics)
library(tictoc)


mixAnalysis <- function(geSimulated, mirnaSimulated, methSimulated, sample_labels, dataname){
  
  sample_labels <- as.factor(sample_labels)
  
  ## 75% of the sample size
  smp_size <- floor(0.75 * nrow(geSimulated))
  
  ## set the seed to make the partition reproductible
  set.seed(600)
  train_ind <- sample(seq_len(nrow(geSimulated)), size = smp_size)

  mrna.train <- geSimulated[train_ind, ]
  mirna.train <- mirnaSimulated[train_ind, ]
  me.train <- methSimulated[train_ind, ]
  Y.train <- sample_labels[train_ind]
  

  data <- list(mirna=mirna.train,
               mrna=mrna.train,
               me=me.train)

  Y <- Y.train
  
  mrna.test <- geSimulated[-train_ind, ]
  mirna.test <- mirnaSimulated[-train_ind, ]
  me.test <- methSimulated[-train_ind, ]
  Y.test <- sample_labels[-train_ind]
  
  data.test <- list(mirna = mirna.test,
                    mrna = mrna.test,
                    me = me.test)

  tic(dataname)
  # creating design matrix and selecting number of components
  design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
  diag(design) = 0
  
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = 3, design = design)
  
  # set.seed(222) 
  perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10, cpus = 24)
  print(perf.diablo$choice.ncomp$WeightedVote)
  
  ER <- "Overall.ER"
  dist <- "mahalanobis.dist"
  
  ncomp = perf.diablo$choice.ncomp$WeightedVote[ER, dist]
  
  test.keepX = list (mirna = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   mrna = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   me = c(5:9, seq(10, 18, 2), seq(20,30,5)))
  set.seed(222)
  tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                               test.keepX = test.keepX, design = design,
                               validation = 'Mfold', folds = 10, nrepeat = 1,
                               dist = dist, cpus = 24)

  list.keepX = tune.TCGA$choice.keepX
  
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design) # final model
  
  
  #performance assessment of the model
  #set.seed(123) # for reproducibility, only when the `cpus' argument is not used
  perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = dist, cpus = 24)
  
  data.test <- list(mirna=mirna.test,
                  mrna=mrna.test,
                  me=me.test)

  predict.diablo = predict(sgccda.res, newdata = data.test)
  toc()
  
  return(list(Y.test, predict.diablo,ncomp))
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

## Small datasets group - Simulated dataset A
load("../../simulated_data/simulatedDatasets/SimulatedCaseA_NN.RData")
mix.ANN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "ANN")

truth <- mix.ANN[[1]]
predict.ANN <- mix.ANN[[2]]
ncomp <- mix.ANN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.ANN[[1]], predicted = predict.ANN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyANN = round(sum(diag(confusion.mat))/length(truth),2)
resANN <- PerformanceAnalysis(confusion.mat)

GroupScoresMix <- data.frame(matrix(,nrow = 24, ncol = 6))
colnames(GroupScoresMix) <- c("Tool", "Dataset", "Group", "Precision", "Recall", "FScore")
GroupScoresMix[,1] <- "mixOmics"
GroupScoresMix[,3] <- rep(c("Group1", "Group2", "Group3"), 8)

DataScoresMix <- data.frame(matrix(,nrow = 8, ncol = 4))
colnames(DataScoresMix) <- c("Tool", "Dataset", "Accuracy", "AvgFScore")
DataScoresMix[,1] <- "mixOmics"

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[1:3,2] <- "A_NN"
GroupScoresMix[1,4:6] <- resANN[[1]]
GroupScoresMix[2,4:6] <- resANN[[2]]
GroupScoresMix[3,4:6] <- resANN[[3]]

DataScoresMix[1,2] <- "A_NN"
DataScoresMix[1,3] <- AccuracyANN
DataScoresMix[1,4] <- resANN[[4]]

## Small noisy datasets group - Simulated dataset A
load("../../simulated_data/simulatedDatasets/SimulatedCaseA_WN.RData")
mix.AWN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "AWN")

truth <- mix.AWN[[1]]
predict.AWN <- mix.AWN[[2]]
ncomp <- mix.AWN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.AWN[[1]], predicted = predict.AWN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyAWN = round(sum(diag(confusion.mat))/length(truth),2)
resAWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[4:6,2] <- "A_WN"
GroupScoresMix[4,4:6] <- resAWN[[1]]
GroupScoresMix[5,4:6] <- resAWN[[2]]
GroupScoresMix[6,4:6] <- resAWN[[3]]

DataScoresMix[2,2] <- "A_WN"
DataScoresMix[2,3] <- AccuracyAWN
DataScoresMix[2,4] <- resAWN[[4]]

## Small datasets group -Simulated dataset B
load("/../../simulated_data/simulatedDatasets/SimulatedCaseB_NN.RData")
mix.BNN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "BNN")

truth <- mix.BNN[[1]]
predict.BNN <- mix.BNN[[2]]
ncomp <- mix.BNN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.BNN[[1]], predicted = predict.BNN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyBNN = round(sum(diag(confusion.mat))/length(truth),2)
resBNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[7:9,2] <- "B_NN"
GroupScoresMix[7,4:6] <- resBNN[[1]]
GroupScoresMix[8,4:6] <- resBNN[[2]]
GroupScoresMix[9,4:6] <- resBNN[[3]]

DataScoresMix[3,2] <- "B_NN"
DataScoresMix[3,3] <- AccuracyBNN
DataScoresMix[3,4] <- resBNN[[4]]

## Small noisy datasets group - Simulated dataset B
load("../../simulated_data/simulatedDatasets/SimulatedCaseB_WN.RData")
mix.BWN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "BWN")

truth <- mix.BWN[[1]]
predict.BWN <- mix.BWN[[2]]
ncomp <- mix.BWN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.BWN[[1]], predicted = predict.BWN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyBWN = round(sum(diag(confusion.mat))/length(truth),2)
resBWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[10:12,2] <- "B_WN"
GroupScoresMix[10,4:6] <- resBWN[[1]]
GroupScoresMix[11,4:6] <- resBWN[[2]]
GroupScoresMix[12,4:6] <- resBWN[[3]]

DataScoresMix[4,2] <- "B_WN"
DataScoresMix[4,3] <- AccuracyBWN
DataScoresMix[4,4] <- resBWN[[4]]

## Small datasets group - Simulated dataset C
load("../../simulated_data/simulatedDatasets/SimulatedCaseC_NN.RData")
mix.CNN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "CNN")

truth <- mix.CNN[[1]]
predict.CNN <- mix.CNN[[2]]
ncomp <- mix.CNN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.CNN[[1]], predicted = predict.CNN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyCNN = round(sum(diag(confusion.mat))/length(truth),2)
resCNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[13:15,2] <- "C_NN"
GroupScoresMix[13,4:6] <- resCNN[[1]]
GroupScoresMix[14,4:6] <- resCNN[[2]]
GroupScoresMix[15,4:6] <- resCNN[[3]]

DataScoresMix[5,2] <- "C_NN"
DataScoresMix[5,3] <- AccuracyCNN
DataScoresMix[5,4] <- resCNN[[4]]

## Small noisy datasets group - Simulated dataset C
load("../../simulated_data/simulatedDatasets/SimulatedCaseC_WN.RData")
mix.CWN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "CWN")

truth <- mix.CWN[[1]]
predict.CWN <- mix.CWN[[2]]
ncomp <- mix.CWN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.CWN[[1]], predicted = predict.CWN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyCWN = round(sum(diag(confusion.mat))/length(truth),2)
resCWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[16:18,2] <- "C_WN"
GroupScoresMix[16,4:6] <- resCWN[[1]]
GroupScoresMix[17,4:6] <- resCWN[[2]]
GroupScoresMix[18,4:6] <- resCWN[[3]]

DataScoresMix[6,2] <- "C_WN"
DataScoresMix[6,3] <- AccuracyCWN
DataScoresMix[6,4] <- resCWN[[4]]

## Small datasets group - Simulated dataset D
load("../../simulated_data/simulatedDatasets/SimulatedCaseD_NN.RData")
mix.DNN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "DNN")

truth <- mix.DNN[[1]]
predict.DNN <- mix.DNN[[2]]
ncomp <- mix.DNN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.DNN[[1]], predicted = predict.DNN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyDNN = round(sum(diag(confusion.mat))/length(truth),2)
resDNN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[19:21,2] <- "D_NN"
GroupScoresMix[19,4:6] <- resDNN[[1]]
GroupScoresMix[20,4:6] <- resDNN[[2]]
GroupScoresMix[21,4:6] <- resDNN[[3]]

DataScoresMix[7,2] <- "D_NN"
DataScoresMix[7,3] <- AccuracyDNN
DataScoresMix[7,4] <- resDNN[[4]]

## Small noisy datasets group - Simulated dataset D
load("../../simulated_data/simulatedDatasets/SimulatedCaseD_WN.RData")
mix.DWN <- mixAnalysis(geSimulated, mirnaSimulated, methSimulated, sample_labels, "DWN")

truth <- mix.DWN[[1]]
predict.DWN <- mix.DWN[[2]]
ncomp <- mix.DWN[[3]]

confusion.mat = get.confusion_matrix(truth = mix.DWN[[1]], predicted = predict.DWN$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat

AccuracyDWN = round(sum(diag(confusion.mat))/length(truth),2)
resDWN <- PerformanceAnalysis(confusion.mat)

## Filling the GroupScores and DataScores table
# GroupScores
GroupScoresMix[22:24,2] <- "D_WN"
GroupScoresMix[22,4:6] <- resDWN[[1]]
GroupScoresMix[23,4:6] <- resDWN[[2]]
GroupScoresMix[24,4:6] <- resDWN[[3]]

DataScoresMix[8,2] <- "D_WN"
DataScoresMix[8,3] <- AccuracyDWN
DataScoresMix[8,4] <- resDWN[[4]]

write.table(DataScoresMix, file = "DataScoresMix.txt", quote = F, sep = "\t")
write.table(GroupScoresMix, file = "GroupScoresMix.txt", quote = F, sep = "\t")
save(mix.ANN, mix.AWN, mix.BNN, mix.BWN, mix.CNN, mix.CWN, mix.DNN, mix.DWN, resANN, resAWN, resBNN, resBWN, resCNN,
     resCWN, resDNN, resDWN, file = "SimulatedSmallDataAnalysis_mixOmics.RData")

