setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/Applications")

source('SideFunctions_Update.R')
source('RidgeFunctions_Update.R')
############ Application ############
load("ClinOmLymphomaApplicationRRChop.Rdata")
rownames(dat) <- NULL
# libraries
library(pROC)
library(dplyr)
library(randomForestSRC)
library(porridge)
library(penalized)
library(magic)
library(rpart); 
library(rpart.plot)
library(treeClust)
library(multiridge)
library(ModTools)
library(caret)
library(xgboost)
ids = which(dat$cohort=="PETAL")
Testdat = dat[ids,]
Traindat = dat[-ids,]
remove(dat,ids)
Ztr = dplyr::select(Traindat, transBCL2,transBCL6,transMYC, IPI, COO)
Zte = dplyr::select(Testdat, transBCL2,transBCL6,transMYC, IPI, COO)
Ytr = Traindat[,145]; Yte = Testdat[,145]
Xtr = Traindat[,1:137]; Xte = Testdat[,1:137]

Xtr = as.matrix(Xtr); Ztr = as.matrix(Ztr)
Xte = as.matrix(Xte); Zte = as.matrix(Zte)
remove(Testdat,Traindat)
set.seed(3)

MSE_ridge <- c(); AUC_ridge <- c()
MSE_FusReg <- c(); AUC_FusReg <- c()
MSE_FullFus <- c(); AUC_FullFus <- c()
MSE_RF <- c(); AUC_RF <- c()
MSE_GB <- c(); AUC_GB <- c()
FusPar <- matrix(NA,nrow=length(Folds),ncol=4)
  
  ######## Model Fitting ######
  
  #1. RF
  DFtrain <- data.frame(Ydf=factor(Ytr),Xdf=cbind(Xtr,Ztr))
  DFtest <- data.frame(Ydf=factor(Yte),Xdf=cbind(Xte,Zte))
  #set.seed(i*2+4*i^2+3)
  RF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=1000, var.used="all.trees",importance=c("none"))
  
  preds_RF <- predict.rfsrc(RF, newdata = DFtest, outcome = "train")$predicted[,2]
  MSE_RF = mean((preds_RF-Yte)^2)
  AUC_RF = roc(Yte,preds_RF)$auc
  remove(DFtrain,DFtest,RF,preds_RF)
  
  #2. Gradient Boosting
  tune_grid <- expand.grid(
    nrounds = c(100,500,1000),
    eta = c(0.01, 0.1, 0.3),
    max_depth = c(3,4),
    gamma = c(0),
    colsample_bytree = 0.75,
    min_child_weight = 20,
    subsample = 1
  )
  
  tune_control <- caret::trainControl(
    method = "cv", # cross-validation
    number = 5, # with n folds 
    verboseIter = FALSE, # no training log
    allowParallel = FALSE #
  )
  
  xgb_tune <- caret::train(
    x = cbind(Xtr,Ztr),
    y = factor(Ytr),
    trControl = tune_control,
    tuneGrid = tune_grid,
    method = "xgbTree",verbosity = 0,verbose=F
  )
  xgb_tune$bestTune
  gbpred <- predict (xgb_tune,cbind(Xte,Zte), type = "prob")[,2]
  MSE_GB = mean((gbpred-Yte)^2)
  AUC_GB = roc(Yte,gbpred)$auc
  remove(xgb_tune,gbpred,tune_control,tune_grid)
  
  # Fit ridge regression
  Lam1=optL2(response = Ytr,penalized = Xtr, unpenalized = cbind(1, Ztr), 
             model = "logistic",trace=F)$lambda
  if (Lam1==Inf){Lam1=10e9}
  RidgeFit = penalized(Ytr,penalized = Xtr,unpenalized = cbind(1, Ztr),lambda2 = Lam1,trace=F, model = "logistic")
  
  
  preds_Ridge = penalized::predict(RidgeFit, penalized=Xte, unpenalized=cbind(1, Zte))
  MSE_ridge = mean((preds_Ridge-Yte)^2);
  AUC_ridge = roc(Yte,preds_Ridge)$auc
  remove(RidgeFit,preds_Ridge)
  
  
  
  # fitting tree
  dat=cbind.data.frame(Ytr,Ztr)
  Treefit <- rpart(Ytr~.,data = dat, control = rpart.control(xval =10, minbucket = 15,cp=0),
              model = T)
  
  prp(Treefit,extra=1)
  idsVar = SplitVariables(Tree = Treefit,Z = Ztr)
  remove(dat)
  
  ####### 4. Fit fully FusedReg Tree ######
  Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,data.frame(Ztr))]
  Nodes <- as.numeric(Nodes)
  
  #sort(unique(Nodes))
  
  names = paste0("N",sort(unique(Nodes)))
  
  Intercepts = model.matrix(~0+factor(Nodes))
  colnames(Intercepts)[1:length(unique(Nodes))] = names
  Un_Train = as.matrix(cbind(Intercepts,Ztr[,idsVar]))
  
  foldsHyp <- CVfolds(Y=Ytr, model = "logistic", balance = T,kfold = 5,nrepeat = 3)
  
  Lam <- optPenaltyGLM.kCVauto(Y = Ytr, X = Xtr, U = Un_Train, 
                               lambdaInit = 10, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(Xtr), nrow = ncol(Xtr)),
                               model ="logistic", folds = foldsHyp, loss = "loglik", 
                               maxIter = 100, implementation = "org")
  
  
  
  Fit_Full <- ridgeGLM(Y = Ytr, U = Un_Train, X = Xtr, 
                       lambda = Lam, lambdaG = 0, Dg = matrix(0, ncol=ncol(Xtr), nrow = ncol(Xtr)), 
                       model="logistic")
  
  Nodes_test <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,data.frame(Zte))]
  Nodes_test <- as.numeric(Nodes_test)
  Intercepts_Test = model.matrix(~0+factor(Nodes_test))
  colnames(Intercepts_Test)[1:length(unique(Nodes))] = names
  Un_Test = as.matrix(cbind(Intercepts_Test,Zte[,idsVar]))
  
  LP = cbind(Un_Test,Xte) %*% Fit_Full
  Ypred_Full =exp(LP)/(1+exp(LP)); Ypred_Full=Ypred_Full[,1]
  MSE_FullFus = mean((Ypred_Full-Yte)^2)
  AUC_FullFus = roc(Yte,Ypred_Full)$auc
  remove(Intercepts,Intercepts_Test,Ypred_Full,LP,names,Nodes,Nodes_test,Fit_Full)
  remove(Un_Test,Un_Train)
  
  
  ####### 5. Fit  FusedReg Tree ######
  Delta = PenMatr2(NumNodes = length(unique(Treefit$where)), p = ncol(Xtr))
  dim(Delta)
  Dat=Dat_Tree(X=Xtr,Z=data.frame(Ztr),Y=Ytr, tree = Treefit, idVars = idsVar)
  Xtr1 <- Dat$X1; Ytr1 <- Dat$Y; Utr1 <- Dat$U
  remove(Dat)
  
  #Estimating penalty parameters lambda and alpha
  begin = proc.time()
  optPenalties <- optPenaltyGLM.kCVauto2(Y = Ytr1, X = Xtr1, U = Utr1, 
                                        lambdaInit = 100, lambdaGinit = 100, Dg=Delta,
                                        model="logistic",
                                        folds=foldsHyp,
                                        loss="loglik", maxIter=100)
  end = proc.time()-begin
  end
  begin = proc.time()
  optPenalties <- porridge::optPenaltyGLM.kCVauto(Y = Ytr1, X = as.matrix(Xtr1), U = Utr1, 
                                         lambdaInit = 100, lambdaGinit = 100, Dg=Delta,
                                         model="logistic",
                                         folds=foldsHyp,
                                         loss="loglik", maxIter=100)
  end1 = proc.time()-begin
  end1
  
  
  
  # Estimate betas
  begin = proc.time()
  Fit = ridgeGLM2(Y = Ytr1, U = Utr1, X = as.matrix(Xtr1), 
                 lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                 model="logistic")
  end = proc.time()-begin
  end
  Fit <- Fit[2:412]
  Fit <- unlist(Fit)
  begin = proc.time()
  Fit1 = porridge::ridgeGLM(Y = Ytr1, U = Utr1, X = as.matrix(Xtr1), 
                  lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                  model="logistic")
  end1 = proc.time()-begin
  end
  end1
  
  all.equal(Fit,Fit1)
  names(Fit) = c(colnames(Utr1),colnames(Xtr1))
  Fit
  # predictions
  Dat_test = Dat_Tree(tree = Treefit, X = Xte, Z = data.frame(Zte), Y = Yte, idVars = idsVar)
  XtestTree = Dat_test$X; YtestTree = Dat_test$Y; UtestTree = Dat_test$U
  remove(Dat_test)
  LP = cbind(UtestTree,XtestTree) %*% as.matrix(Fit)
  LP <- LP[,1]
  Ypred =exp(LP)/(1+exp(LP))
  MSE_FusReg = mean((Ypred-YtestTree)^2)
  AUC_FusReg = roc(YtestTree,Ypred)$auc
  FusPar[i,] = c(optPenalties, Lam,Lam1)
  remove(XtestTree,YtestTree,UtestTree,Fit,foldsHyp,Delta,optPenalties,Lam,Lam1,idsVar,Ypred,LP)
  remove(Xtr,Xte,Ztr,Zte,Ytr,Yte,Xtr1,Utr1,Ytr1)
}
res <- cbind.data.frame(AUC_FullFus,AUC_FusReg,AUC_GB,AUC_RF,AUC_ridge,
                        MSE_FullFus,MSE_FusReg,MSE_GB,MSE_RF,MSE_ridge,FusPar)  
save(res, file = "LymphApp_Rchop.Rdata")
round(colMeans(res),2)
