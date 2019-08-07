## GBM data analysis: Analysis with mixOmics
rm(list = ls())
library(mixOmics)
library(tictoc)

gbm.dat.me = readRDS("GBM_Methylation_Top25Filtered.rd")
gbm.dat.rna = readRDS("GBM_RNA_Top25Filtered.rd")
gbm.dat.protein = readRDS("GBM_Protein_Top25Filtered.rd")
data.labels = readRDS("GBM_DataLabels.rd")

data.labels$Label <- as.factor(data.labels$Label)
table(data.labels$Label)

## Creating train and test dataset
## 75% of the sample size
smp_size <- floor(0.75 * nrow(gbm.dat.protein))

## set the seed to make the partition reproductible
set.seed(300)
train_ind <- sample(seq_len(nrow(gbm.dat.protein)), size = smp_size)

protein.train <- gbm.dat.protein[train_ind, ]
mrna.train <- gbm.dat.rna[train_ind, ]
me.train <- gbm.dat.me[train_ind, ]
Y.train <- data.labels[train_ind,]

protein.test <- gbm.dat.protein[-train_ind, ]
mrna.test <- gbm.dat.rna[-train_ind, ]
me.test <- gbm.dat.me[-train_ind, ]
Y.test <- data.labels[-train_ind, ]

data <- list(protein=protein.train,
             mrna=mrna.train,
             me=me.train)
lapply(data, dim)

Y <- Y.train$Label
summary(Y)

## running mixOmics
tic("GBM-mix")
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

sgccda.res = block.splsda(X = data, Y = Y, ncomp=3, design = design)

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


selectVar(sgccda.res, block = 'protein', comp = 1)$protein$name 
selectVar(sgccda.res, block = 'mrna', comp = 1)$mrna$name 
selectVar(sgccda.res, block = 'me', comp = 1)$me$name 


perf.diablo = perf(sgccda.res, validation = 'Mfold', M = 10, nrepeat = 10, dist = 'mahalanobis.dist', cpus = 24)
perf.diablo$MajorityVote.error.rate

data.test <- list(protein=protein.test,
                  mrna=mrna.test,
                  me=me.test)

predict.diablo = predict(sgccda.res, newdata = data.test)
toc()

confusion.mat = get.confusion_matrix(truth = Y.test$Label, predicted = predict.diablo$WeightedVote$mahalanobis.dist[,ncomp])
confusion.mat
get.BER(confusion.mat)

Accuracy = round(sum(diag(confusion.mat))/length(Y.test$Label),4)

## Classical
Pg1 <- round(confusion.mat[1,1]/sum(confusion.mat[,1]),2) # Precision
Rg1 <- round(confusion.mat[1,1]/sum(confusion.mat[1,]),2) # Recall
Fg1 <- round(2*((Pg1*Rg1)/(Pg1+Rg1)),2) # F1-score
Fg1[is.nan(Fg1)] <- 0

g1 <- c(Pg1, Rg1, Fg1)

# Mesenchymal
Pg2 <- round(confusion.mat[2,2]/sum(confusion.mat[,2]),2) # Precision
Rg2 <- round(confusion.mat[2,2]/sum(confusion.mat[2,]),2) # Recall
Fg2 <- round(2*((Pg2*Rg2)/(Pg2+Rg2)),2) # F1-score
Fg2[is.nan(Fg2)] <- 0

g2 <- c(Pg2, Rg2, Fg2)

# Neural
Pg3 <- round(confusion.mat[3,3]/sum(confusion.mat[,3]),2) # Precision
Rg3 <- round(confusion.mat[3,3]/sum(confusion.mat[3,]),2) # Recall
Fg3 <- round(2*((Pg3*Rg3)/(Pg3+Rg3)),2) # F1-score
Fg3[is.nan(Fg3)] <- 0

g3 <- c(Pg3, Rg3, Fg3)

# Proneural
Pg4 <- round(confusion.mat[4,4]/sum(confusion.mat[,4]),2) # Precision
Rg4 <- round(confusion.mat[4,4]/sum(confusion.mat[4,]),2) # Recall
Fg4 <- round(2*((Pg4*Rg4)/(Pg4+Rg4)),2) # F1-score
Fg4[is.nan(Fg4)] <- 0

g4 <- c(Pg4, Rg4, Fg4)


AvgF <- round((Fg1+Fg2+Fg3+Fg4)/4,2)
save.image(file= "GBM_mixOmicsAnalysis_withPerformanceScore_Top25Filtered.RData")

