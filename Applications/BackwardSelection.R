######### Backward Selection #########

# check whether omics add to leaf nodes
setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Applications")
library(survival)
library(survC1)
library(survminer)
library(rpart); 
library(rpart.plot)
library(treeClust)
library(multiridge)
library(penalized)
library(glmnet)
library(globaltest)
library(randomForestSRC)
library(blockForest)
library(gbm)
library(SurvMetrics)
library(survivalROC)
library(stringr)

load("FinalFitCRC.Rdata")
remove(U1,X2,Fit,Y)

load("CRCdata.Rdata")
set.seed(48)
ids=sample(1:nrow(Clinical),size = 0.2*nrow(Clinical),replace = F)
X <-Omics[-ids,]; Xtest <- Omics[ids,]
Y <- Response[-ids,]; Ytest <- Response[ids,]
Z<-Clinical[-ids,]; Ztest=Clinical[ids,]
remove(ids,Response,Clinical,Omics)
gc()


rpart.plot(Treefit, 
           type = 5,
           extra = 1, 
           box.palette = "Pu",
           branch.lty = 8, 
           shadow.col= 0, 
           nn = TRUE,
           cex = 0.6)
p_clin = ncol(model.matrix(~.,Z)[,-1])
p = ncol(X)


Lin = T
RelTime = 5
tau = 8

set.seed(5)
foldsHyp = CVfoldsTree(Y = Y, Tree = Treefit, Z = Z, model = "cox", kfold = 5, nrepeat = 3)
Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
Nodes <- as.numeric(Nodes)
NumNod <- length(unique(Nodes))

Nodenames <- paste0("N",sort(unique(Nodes)))
pvals = sapply(sort(unique(Nodes)), function(x) globaltest::p.value(globaltest::gt(
  Surv(Y[Nodes==x,1],Y[Nodes==x,2]), 
  alternative = data.frame(matrix(X[Nodes == x, ],nrow =length(which(Nodes==x)),ncol=ncol(X))),
  model = "cox")))
names(pvals) <- Nodenames
pvals
pvals[is.na(pvals)]<-1
print("Estimated p values in the nodes are:")
print(pvals)
pvals<-sort(pvals,decreasing = T)
print("Remove omics effects in nodes in order:")
print(paste(names(pvals),collapse = " -> "))
Dat = Dat_Tree(Tree = Treefit, X = X, Z = Z,  LinVars = Lin)
X1 = Dat$Omics; U1 = Dat$Clinical
colSums(U1)
Dat_tree_test = Dat_Tree(Tree=Treefit, X=Xtest, Z=Ztest, LinVars = Lin)
Ute_tree = Dat_tree_test$Clinical; Xte_tree = Dat_tree_test$Omics

remove(Dat,Dat_tree_test,Yte_tree)
Concordances1 <- c()
AUCs <- c()
for (j in 1:length(pvals)) {
  j = 3
  EmpNodes = names(pvals)[1:j]
  #names(Fits)[i] <- paste(names(pvals)[1:i],collapse = ",")
  print(paste("Fit FusedTree without omics effects in", paste(names(pvals)[1:j],collapse = ", ")))
  ids = lapply(EmpNodes, function (x) grepl(x,str_sub(colnames(X1), start = -3)))
  ids1 = which(Reduce("|",ids))
  colnames(X1)[1:15]
  X2 <- X1[,-ids1]
  Xte <- Xte_tree[,-ids1]
  dim(X2)
  remove(ids,ids1)
  NodNum <- ncol(X2)/p
  print(NodNum)
  
  if (NodNum == 0){
    ClinFit <- coxph(Y ~.,data = data.frame(U1[,-1]))
    summary(ClinFit)
    LP <- predict(ClinFit,newdata = Z,type = "lp")
    LpPred <- predict(ClinFit,newdata = data.frame(Ute_tree[,-1]) ,type = "lp")
    #ConcClin <- concordance(Ytest ~ LpPred)$concordance
    Concordances1[j] <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LpPred), tau = tau,  nofit=T)$Dhat
    AUCs[j] <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LpPred, predict.time = RelTime)$AUC
    
    remove(Fit,ids)
    
  }
  if (NodNum == 1){
    dim(X2)
    
    Lambda <- optPenaltyGLM.kCVauto2(Y = Y1, X = X2, U = U1, 
                                     lambdaInit = Lam1, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X2), nrow = ncol(X2)),
                                     folds = foldsHyp,
                                     model = "surv", loss="loglik")
    print(Lambda)
    Fit <-ridgeGLM2(Y = Y1, U = U1, X = X2, 
                    lambda = Lambda, lambdaG = 0, Dg = matrix(0, ncol=ncol(X2), nrow = ncol(X2)), 
                    model="surv")
    names(Fit)<-c(colnames(U1),colnames(X2))
    LP = -as.numeric(cbind(Ute_tree,Xte) %*% Fit)
    Concordances1[j] <- Est.Cval(cbind(Ytest[,1],Ytest[,2],-LP), tau=tau, nofit=T)$Dhat
    AUCs[j] <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = -LP, predict.time = RelTime)$AUC
    remove(Fit,LP,Preds)
    
  }
  if (NodNum > 1){
    
    Delta = .PenMatr(NumNodes = NodNum, p = p)
    dim(Delta)
    dim(X2)
    gc()
    if (ncol(Delta) != ncol(X2)){stop("number of columns of penalty matrix does not equal number of columns of design matrix X")}
    #if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
    #if(alphaInit <= 0){stop("initial penalty alphaInit should be larger than zero")}
    optPenalties <- .optPenaltyGLM.kCVauto2(Y = Y, X = X2, U = U1, 
                                           lambdaInit = 1500, lambdaGinit = 100,
                                           folds = foldsHyp,
                                           Dg = Delta,model = "cox", loss="loglik")
    
    optPenalties <- .optPenaltyGLM.kCVauto2(Y = Y, X = X2, U = U1, 
                                           lambdaInit = 1500, lambdaGinit = 0,
                                           folds = foldsHyp,
                                           model = "cox",loss="loglik")
    
    print(optPenalties)
    Fit = .ridgeGLM2(Y = Y, U = U1, X = X2, 
                    lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                    model="cox")
    Fit = .ridgeGLM2(Y = Y, U = U1, X = X2, 
                    lambda = optPenalties[1], lambdaG = 0,
                    model= "cox", verbose = TRUE)
    Fit <- Fit$estimates
    names(Fit)<-c(colnames(U1),colnames(X2))
    LP = as.numeric(cbind(Ute_tree,Xte) %*% Fit)
    Concordances1[j] <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LP), tau=tau, nofit=T)$Dhat
    AUCs[j] <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP, predict.time = RelTime)$AUC
    remove(Fit)
  }
  remove(X2)
}
Concordances1
AUCs

is.numeric(Y[-foldsHyp[[2]], 1])
