## hepatocellular carcinoma dataset 2 labels: Analysis with mixOmics
rm(list = ls())
library(mixOmics)
library(tictoc)

load("/home/n9870130/Benchmark_Integrative_Tools/HepatocellularCarcinoma/HPC_filtered_data.Rdata")

# Log tranforming the data as they are in counts
mirna2 <- log(mirna+1)
mrna2 <- log(mrna+1)
rm(mirna, mrna)

data.labels$Labels <- as.factor(data.labels$Labels)
table(data.labels$Labels)

## Creating train and test dataset
## 75% of the sample size
smp_size <- floor(0.75 * nrow(mirna2))

## set the seed to make the partition reproductible
set.seed(300)
train_ind <- sample(seq_len(nrow(mirna2)), size = smp_size)

mirna.train <- mirna2[train_ind, ]
mrna.train <- mrna2[train_ind, ]
me.train <- me[train_ind, ]
Y.train <- data.labels[train_ind,]

mirna.test <- mirna2[-train_ind, ]
mrna.test <- mrna2[-train_ind, ]
me.test <- me[-train_ind, ]
Y.test <- data.labels[-train_ind, ]

data <- list(mirna=mirna.train,
             mrna=mrna.train,
             me=me.train)
lapply(data, dim)

Y <- Y.train$Labels
summary(Y)

## running mixOmics
tic("HPC3Labels-mix")
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

sgccda.res = block.splsda(X = data, Y = Y, ncomp=2, design = design)

perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 10, nrepeat = 10, cpus = 24)

perf.diablo$choice.ncomp$WeightedVote
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.ER", "mahalanobis.dist"]

test.keepX = list (mirna = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   mrna = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   me = c(5:9, seq(10, 18, 2), seq(20,30,5)))


tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                               test.keepX = test.keepX, design = design,
                               validation = 'Mfold', folds = 10, nrepeat = 1,
                               cpus = 24, dist = "mahalanobis.dist")
list.keepX = tune.TCGA$choice.keepX

sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design) # final model


selectVar(sgccda.res, block = 'mirna', comp = 1)$mirna$name 
selectVar(sgccda.res, block = 'mrna', comp = 1)$mrna$name 
selectVar(sgccda.res, block = 'me', comp = 1)$me$name 


perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'mahalanobis.dist', cpus = 24)
perf.diablo$MajorityVote.error.rate

data.test <- list(mirna=mirna.test,
                  mrna=mrna.test,
                  me=me.test)

predict.diablo = predict(sgccda.res, newdata = data.test)
toc()

confusion.mat = get.confusion_matrix(truth = Y.test$Labels, predicted = predict.diablo$WeightedVote$mahalanobis.dist)
confusion.mat
get.BER(confusion.mat)

Accuracy = round(sum(diag(confusion.mat))/length(Y.test$Labels),4)

## Normal
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0
  
g1 <- c(Pg1, Rg1, Fg1)

# PVTT
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

g2 <- c(Pg2, Rg2, Fg2)

# tumour
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0
  
g3 <- c(Pg3, Rg3, Fg3)

AvgF <- round((Fg1+Fg2+Fg3)/3,2)

save.image(file= "mixOmics_HPC3Labels.RData")
