setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Applications")

library(survival)
library(survC1)
library(survminer)
library(rpart) 
library(rpart.plot)
library(treeClust)
library(multiridge)
library(globaltest)
library(SurvMetrics)
library(survivalROC)
library(stringr)
#remotes::install_github("JeroenGoedhart/fusedTree")
library(fusedTree)

load("CRCdata.Rdata")

set.seed(48)
ids=sample(1:nrow(Clinical),size = 0.2*nrow(Clinical),replace = F)
Y<- Surv(time = Response[,1], event = Response[,2])
X<-Omics[-ids,]; Xtest<-Omics[ids,]
Y<-Response[-ids,]; Ytest<-Response[ids,]
Z<-Clinical[-ids,]; Ztest=Clinical[ids,]
remove(ids,Response,Clinical,Omics)
gc()

p_clin = ncol(model.matrix(~.,Z)[,-1])
p = ncol(X)
Lin = T
RelTime = 5
tau = 8

load("FinalFitCRC.Rdata")
remove(U1,X2); gc()

rpart.plot(Treefit, # middle graph
           type=5,
           extra=1, 
           box.palette="Pu",
           branch.lty=8, 
           shadow.col="gray", 
           nn=TRUE,
           cex = 0.6)

Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
Nodes <- as.numeric(Nodes)
NumNod <- length(unique(Nodes))
set.seed(5)
foldsHyp <- CVfoldsTree(Y=Y,Tree = Treefit,Z=Z,model = "cox", 
                        kfold = 5, nrepeat = 3)

optLam1 <- PenOpt(Tree = Treefit,Y = Y, X = X, Z = Z, model = "cox",
                  lambdaInit = 1000, alphaInit = 0,
                  folds = foldsHyp,LinVars = Lin, maxIter = 30)
optLam1

Fit <- fusedTree(Tree=Treefit,Y = Y, X = X, Z = Z,
                  model = "cox",
                  lambda = optLam1[1],
                  alpha =  0,
                  LinVars = Lin)
Preds <- predict(Fit,newX = Xtest,newZ = Ztest, newY = Ytest)
Preds <- Preds$LinPred[,2]

ConcFusTree <- Est.Cval(cbind(Ytest[,1], Ytest[,2], Preds), 
                     tau = tau, nofit = TRUE)$Dhat
AUCFusTree <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], 
                            marker = Preds, predict.time = RelTime)$AUC

effects <- Fit$Effects
eff_clin <- effects[1:7]
eff_omics <- effects[-(1:7)]

Nodenames <- paste0("N",sort(unique(Nodes)))
pvals <- sapply(sort(unique(Nodes)), function(x) globaltest::p.value(globaltest::gt(
  Surv(Y[Nodes == x, 1], Y[Nodes == x, 2]), 
  alternative = data.frame(matrix(X[Nodes == x, ], nrow = length(which(Nodes == x)),ncol = ncol(X))),
  model = "cox")))
names(pvals) <- Nodenames
pvals
pvals[is.na(pvals)] <- 1
print("Estimated p values in the nodes are:")
print(pvals)
print("Remove omics effects in nodes in order:")
pvals <- sort(pvals, decreasing = TRUE)
print(paste(names(pvals),collapse = " -> "))
#Dat <- Dat_Tree(Tree = Treefit, X = X, Z = Z,  LinVars = Lin)
#U <- Dat$Clinical; X1 <- Dat$Omics
Dat_tree_test <- Dat_Tree(Tree=Treefit, X = Xtest, Z = Ztest, LinVars = Lin)
U_test <- Dat_tree_test$Clinical; X1_test <- Dat_tree_test$Omics
remove(Dat_tree_test)
gc()
Concordances <- c()
AUCs <- c()
for (j in 1:length(pvals)) {
  
  EmpNodes <- names(pvals)[1:j]
  #names(Fits)[i] <- paste(names(pvals)[1:i],collapse = ",")
  print(paste("Fit FusedTree without omics effects in", paste(names(pvals)[1:j],collapse = ", ")))
  ids <- lapply(EmpNodes, function (x) grepl(x,str_sub(colnames(X1), start = -3)))
  ids1 <- which(Reduce("|",ids))
  X2 <- X1[,-ids1]
  Xte <- X1_test[,-ids1]
  print(dim(Xte))
  LP <- as.numeric(cbind(U_test,Xte) %*% c(eff_clin,eff_omics[-ids1]))
  Concordances[j] <- Est.Cval(cbind(Ytest[,1], Ytest[,2], LP), 
                               tau = tau, nofit = TRUE)$Dhat
  AUCs[j] <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP, 
                           predict.time = RelTime)$AUC
  remove(ids, ids1, X2, Xte, LP, EmpNodes)
}
AUCs
