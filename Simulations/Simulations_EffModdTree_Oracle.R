setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/Simulations")

source('SideFunctions_Update.R')
source('RidgeFunctions_Update.R')
#libraries

library(porridge)
library(Matrix)
library(magic)
library(rpart); 
library(rpart.plot)
library(treeClust)
library(multiridge)
library(ModTools)
library(mvtnorm)
library(globaltest)
colMeans(res)
##### simulations #####

Nsim = 100
p = 500
p_Clin = 5
N = 300
Ntest = 5000
load("C:/Users/VNOB-0732/Desktop/R files/Learn2Evaluate/Learn2Evaluate/Simulation Results/Data/correlationmatrix.Rdata")
CorrX <- CorrX[1:p,1:p]


set.seed(23)
Xtest = mvtnorm::rmvnorm(Ntest, sigma = CorrX)
colnames(Xtest) = paste0("x",seq(1,p))
Ztest = matrix(rnorm(Ntest*p_Clin,0,1),nrow = Ntest, ncol = p_Clin)
colnames(Ztest) = paste0("z",seq(1,p_Clin))
Ztest = data.frame(Ztest)
#Var(Xtest[,1:25]%*%betas[1:25])

g <- function (z,x,betas){
    1*(z[,1] < 0)*(z[,2] < 0)*(-10 + x[,1:p/2] %*% betas[1:p/2]*4) + 
    1*(z[,1] < 0)*(z[,2] >= 0)*(-5 + x[,1:p/2]%*% betas[1:p/2]*2) +
    1*(z[,1] >= 0)*(z[,4] < 0)*(5 + x[,1:p/2]%*% betas[1:p/2]*0.5) +
    1*(z[,1] >= 0)*(z[,4] >= 0)*(10 + x[,1:p/2]%*% betas[1:p/2]*0.1) +
    x[,(p/2+1):p]%*% betas[(p/2+1):p] 
}



treef <- function(z){
  nodes <- ifelse((z[,1] < 0 & z[,2] < 0),1,0) +
           ifelse((z[,1] < 0 & z[,2] >= 0),2,0) +
           ifelse((z[,1] >= 0 & z[,4] < 0),3,0) +
           ifelse((z[,1] >= 0 & z[,4] >= 0),4,0)
  
  return(nodes)
}
Nodes_or_Test = treef(Ztest)

unique(Nodes_or_Test)
summary(factor(Nodes_or_Test))

set.seed(1)
betas<- matrix(rexp(p,p/10))
Ytest = g(z = Ztest, x = Xtest, betas=betas)[,1] + rnorm(Ntest,0,1)
var(Ytest)


Nsim =100

MSE_FusReg <- c()
MSE_FullFus <- c()
MSE_Oracle <- c()

FusPar <- matrix(NA,nrow=Nsim,ncol=5)

colnames(FusPar) <- c("Lambda_or","alpha_or","Lambda","alpha","Lam Full Fus")
VarsBeta <- matrix(NA,nrow=Nsim,ncol=2)

for (i in 1:Nsim) {
  
  print(paste("simulation",i, sep = " "))
  ## simulating clinical covariates
  
  set.seed(i^3+52)
  Z = matrix(rnorm(N*p_Clin,0,1),nrow = N, ncol = p_Clin) 
  colnames(Z) = paste0("z",seq(1,p_Clin))
  Z <- data.frame(Z)
  X = mvtnorm::rmvnorm(N, sigma = CorrX)
  colnames(X) = paste0("x",seq(1,p))
  
  ## simulating response
  #set.seed(i*3+4*i^2+3)
  Y = g(z = Z, x = X, betas=betas)[,1] + rnorm(N,0,1)
  var(Y)
  ####### 3. Fit Ridge regression (Clinical unp/penalized) ######
  folds <- CVfolds(Y=Y, model = "linear", balance = F,kfold = 5,nrepeat = 1)
  
  
  ### oracle tree model
  Nodes_oracle <- treef(Z)
  summary(factor(Nodes_oracle))
  Delta_or = PenMatr2(NumNodes = length(unique(Nodes_oracle)), p = p)
  
  Dat=Dat_Tree(X=X,Z=Z,Y=Y, nodes = Nodes_oracle, idVars = 0)
  X_or <- Dat$X1; Y_or <- Dat$Y; U_or <- Dat$U
  
  remove(Dat)
  optPenalties_or <- optPenaltyGLM.kCVauto2(Y = Y_or, X = X_or, U = U_or, 
                                         lambdaInit = 100, lambdaGinit = 1, Dg=Delta_or,
                                         model="linear",
                                         folds=folds,
                                         loss="sos", maxIter=100)
  
  Fit_or <- ridgeGLM(Y = Y_or, U = U_or, X = X_or, 
                     lambda = optPenalties_or[1], lambdaG = optPenalties_or[2], Dg = Delta_or, 
                     model="linear")
  
  
  
  Dat_test = Dat_Tree(nodes = Nodes_or_Test, X = Xtest, Z = Ztest, Y = Ytest, idVars = 0)
  XtestTree_or = Dat_test$X1; YtestTree_or = Dat_test$Y; UtestTree_or = Dat_test$U
  remove(Dat_test)
  Ypred = cbind(UtestTree_or,XtestTree_or) %*% as.matrix(Fit_or)
  Ypred <- Ypred[,1]
  MSE_Oracle[i] =mean((Ypred-YtestTree_or)^2)
  remove(XtestTree_or,YtestTree_or,UtestTree_or,Ypred,Fit_or,X_or,Y_or,U_or,Nodes_oracle,Delta_or)
  
  # fitting tree
  dat=cbind.data.frame(Y,Z)
  rp <- rpart(Y~.,data = dat, control = rpart.control(xval =10, minbucket = 30),
              model = T)
  cp = rp$cptable[,1][which.min(rp$cptable[,4])]
  Treefit <- prune(rp, cp = cp)
  prp(Treefit,extra=1)
  #idsVar = SplitVariables(Tree = Treefit, Z= data.frame(Z))
  remove(dat,rp)
  
  ####### 4. Fit fully FusedReg Tree ######
  Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
  Nodes <- as.numeric(Nodes)
  Intercepts = model.matrix(~0+factor(Nodes))
  #sort(unique(Nodes))
  
  names = paste0("N",sort(unique(Nodes)))
  if (length(unique(Nodes))<2){
    Intercepts <- Z
    Intercepts <- model.matrix(~.,Intercepts)
  } else {
    Intercepts = model.matrix(~0+factor(Nodes))
    colnames(Intercepts)[1:length(unique(Nodes))] = names
  }
 
  
  #Un_Train = as.matrix(cbind(Intercepts,Z[,idsVar]))
  
  
  Lam <- porridge::optPenaltyGLM.kCVauto(Y = Y, X = X, U = Intercepts, 
                                         lambdaInit = 100, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)),
                                         model ="linear", folds = folds, loss = "sos", 
                                         maxIter = 100)
  
  
  
  Fit_Full <- ridgeGLM(Y = Y, U = Intercepts, X = X, 
                       lambda = Lam, lambdaG = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)), 
                       model="linear")
  Nodes_test <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Ztest)]
  Nodes_test <- as.numeric(Nodes_test)
  
  if (length(unique(Nodes))<2){
    Intercepts_Test <- Ztest
    Intercepts_Test <- model.matrix(~.,Intercepts_Test)
  } else {
    Intercepts_Test = model.matrix(~0+factor(Nodes_test))
    colnames(Intercepts_Test)[1:length(unique(Nodes))] = names
  }
  Ypred_Full = cbind(Intercepts_Test,Xtest) %*% Fit_Full
  MSE_FullFus[i] = mean((Ypred_Full[,1]-Ytest)^2); MSE_FullFus/var(Ytest)
  var(Ytest)
  remove(Intercepts,Intercepts_Test,Ypred_Full,Fit_Full)
  #remove(Un_Test,Un_Train)
  
  
  ####### 5. Fit  FusedReg Tree ######
  if (length(unique(Nodes))<2){
    MSE_FusReg[i] <- MSE_FullFus[i]
    FusPar[i,] = c(optPenalties_or,Lam,10e10, Lam)
    VarsBeta[i,] = rep(NA,2)
  } else {
    Delta = .PenMatr(NumNodes = length(unique(Nodes)), p = p)
    
    Dat=Dat_Tree(X=X,Z=Z,Y=Y, tree = Treefit,LinVars = F, globaltest = F, model = "linear", p_select = 0.05)
    
    
    #Estimating penalty parameters lambda and alpha
    
    
    optPenalties <- optPenaltyGLM.kCVauto2(Y = Y, X = X3, U = U, 
                                           lambdaInit = 100, lambdaGinit = 1, Dg=as.matrix(Delta),
                                           model="linear",
                                           folds=folds,
                                           loss="sos", maxIter=100)
    
    #begin = proc.time()
    #optPenalties1 <- porridge::optPenaltyGLM.kCVauto(Y = Y, X = X, U = U, 
    #                                      lambdaInit = 10, lambdaGinit = 10, Dg=Delta,
    #                                     model="linear",
    #                                    folds=folds,
    #                                   loss="sos", maxIter=100,
    #                                  implementation="alt")
    #end1 = proc.time()-begin
    #end1
    
    
    # Estimate betas
    Fit = ridgeGLM(Y = Y, U = U, X = X, 
                   lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                   model="linear")
    names(Fit) = c(names,colnames(X))
    # predictions
    Dat_test = Dat_Tree(nodes = Nodes_test, X = Xtest, Z = Ztest, Y = Ytest, idVars = 0)
    XtestTree = Dat_test$X1; YtestTree = Dat_test$Y; UtestTree = Dat_test$U
    remove(Dat_test)
    Ypred = cbind(UtestTree,XtestTree) %*% as.matrix(Fit)
    Ypred <- Ypred[,1]
    MSE_FusReg[i] =mean((Ypred-YtestTree)^2)
    Fit1 = Fit[-(1:ncol(U))]
    names(Fit1) <- colnames(XtestTree)
    Fit1 <- split(Fit1,ceiling(seq_along(Fit1)/ncol(U)))
    Vars = unlist(lapply(Fit1, var))
    FusPar[i,] = c(optPenalties_or,optPenalties, Lam)
    VarsBeta[i,] = c(mean(Vars[1:250]),mean(Vars[251:500]))
    remove(XtestTree,YtestTree,UtestTree,Fit,folds,Delta,optPenalties,optPenalties_or,Lam,Ypred,Fit1,Vars,names)
    remove(Nodes,Nodes_test)
    remove(X,Y,U,Z,Treefit)
    gc()
  }
}

var(Ytest)
View(FusPar)
apply(FusPar, 2, median)
mean(MSE_FusReg);mean(MSE_FullFus)

mean(MSE_Oracle)
plot(FusPar[,2],FusPar[,4])
View(FusPar)
var(Ytest)
results = cbind.data.frame(MSE_FullFus,MSE_FusReg,MSE_GB,MSE_RF,MSE_ridge)#,FusPar,VarsBeta)
round(colMeans(results),1)
round(apply(results,MARGIN = 2, median),3)

nm = paste(p,p_Clin,N,Nsim,"EffModd.Rdata",sep = "_")
save(results,file = nm)




var(Ytest)

mean(MSE_FullFus)
mean(MSE_FusReg)
mean(MSE)
results1 <- results
colMeans(results)
colMeans(results1)
