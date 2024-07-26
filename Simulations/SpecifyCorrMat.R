library(corpcor)
library(glmnet)
library(pROC)
library(randomForest)
library(mvtnorm)

setwd("C:/Users/VNOB-0732/Desktop/R files/Learning Curves/Simulation/R scripts/Simulation-Study-Learning-Curves")
#### defining the correlation matrix, test set and beta's.


### The data simulation set-up ###
Nsim = 1000; #number of simulations, for each simulation a confidence bound is returned
n_train = 100; #number of training samples, will be generated Nsim times  
p = 2000; #number of covariates
n_test = 2.5e4 #test data, will be generated only once



# define the correlation  for covariates. This is done by estimating the correlation matrix
#of a subset (such that p =2000) of a omics dataset.
#then cholesky decomposition is used to generate samples with the desired correlation structure

#preproccesing the omics dataset
load("C:/Users/VNOB-0732/Desktop/R files/Learn2Evaluate/Learn2Evaluate/Data/Bloodplatelet_RNAseq.Rdata")
grnames <- levels(group)
compareGroups = c("NSCLC","HC")
g1 = compareGroups[1]; g2 = compareGroups[2]
id1 = which(group==g1) 
id2 = which(group==g2)
datSqrt2 = dataSqrt[,c(id1,id2)]
sds = apply(datSqrt2,1,sd)
id.del = which(sds==0)
datSqrt2 = datSqrt2[-id.del,]

#standardize data
datStd = t(apply(datSqrt2,1,function(x){(x-mean(x))/sd(x)}))
datStd = t(datStd)
remove(dataSqrt, datSqrt2, g1, g2, group, id1, id2,id.del,sds,compareGroups,grnames)
colnames(datStd) <- NULL

set.seed(1)
X_covariate <- datStd[1:115,c(sample(1:ncol(datStd),p))] # use first p covariates to estimate covariance matrix
dim(X_covariate)
CovX <- cov.shrink(X_covariate) ##estimate covariance matrix by using shrinkage method of Schafer et al (2005)


CorrX[1:10,1:10]
CorrX <-as.matrix.data.frame(CovX)
CorrX <- as.matrix(CorrX)
names(CorrX) <-NULL

isSymmetric.matrix(CorrX) #check if matrix is symmetric
is.positive.definite(CorrX) #check if matrix is positive definite

#saving correlation structure
filenm <- paste("correlationmatrix",".Rdata", sep="")
save(CorrX,file = filenm)