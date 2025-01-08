setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/FusedTree_Paper/Applications")
load("500_CRC.Rdata")
load("5000_CRC.Rdata")
load("21292_CRC.Rdata")
dim(Omics)
gc()
RelTime = 5
tau = 8
# 0. Clinical

library(survival)
library(survC1)
library(multiridge)
library(penalized)
library(SurvMetrics)
library(survivalROC)


# standardize data
Omics <- scale(Omics)
#set.seed(4)
summary(Clinical)
# split data in training and test set
set.seed(48)
#set.seed(163)
#set.seed(4889)

ids=sample(1:nrow(Clinical),size = 0.2*nrow(Clinical),replace = F)
Y<- Surv(time = Response[,1], event = Response[,2])
X<-Omics[-ids,]; Xtest<-Omics[ids,]
Y<-Response[-ids,]; Ytest<-Response[ids,]
Z<-Clinical[-ids,]; Ztest=Clinical[ids,]
remove(ids,Response,Clinical,Omics)



gc()
LP_Tot <- c()
AUCs <- c()
Concs <- c()

for (i in 1:4) { # 4 is number of stage categories
  
  Ys <- split(Y, f = Z$Stage)[[i]]
  Xs <- as.matrix(split(data.frame(X), f = Z$Stage)[[i]])
  Zs <- split(Z[,-4], f = Z$Stage)[[i]]
  Ys_Te <- split(Ytest, f = Ztest$Stage)[[i]]
  Xs_Te <- as.matrix(split(data.frame(Xtest), f = Ztest$Stage)[[i]])
  Zs_Te <- split(Ztest[,-4], f = Ztest$Stage)[[i]]
  Zs <- model.matrix(~.,Zs)[,-1] # create design matrix
  Zs_Te <- model.matrix(~.,Zs_Te)[,-1]
  
  if (i == 1){
    # perfect seperation so fit penalized reg with very small penalty
    RidgeFit <- penalized::penalized(response = Ys, penalized = Zs,
                                     model = "cox", lambda2=0.0001)
    coefsPenalized <- penalized::coefficients(RidgeFit)
    LP <- (Zs_Te %*% coefsPenalized)[,1]
    LP_Tot <- append(LP_Tot,LP)
    Concs[i] <- Est.Cval(cbind(Ys_Te[,1],Ys_Te[,2],LP), tau = tau,  nofit=T)$Dhat
    AUCs[i] <- survivalROC.C(Stime = Ys_Te[,1], status = Ys_Te[,2], marker = LP, predict.time = RelTime)$AUC
    remove(LP,RidgeFit,coefsPenalized,Ys,Zs,Xs,Xs_Te,Ys_Te,Zs_Te)
    gc()
  } else {
  Lam <- penalized::optL2(response = Ys, penalized = Xs, unpenalized = Zs,
                          model = "cox", fold = 10)
  RidgeFit <- penalized::penalized(response = Ys, penalized = Xs, unpenalized = Zs,
                                   model = "cox", lambda2=Lam$lambda)
  coefsPenalized <- penalized::coefficients(RidgeFit)
  LP <- (cbind(Zs_Te,Xs_Te) %*% coefsPenalized)[,1]
  LP_Tot <- append(LP_Tot,LP)
  Concs[i] <- Est.Cval(cbind(Ys_Te[,1],Ys_Te[,2],LP), tau = tau,  nofit=T)$Dhat
  AUCs[i] <- survivalROC.C(Stime = Ys_Te[,1], status = Ys_Te[,2], marker = LP, predict.time = RelTime)$AUC
  remove(LP,Lam,RidgeFit,coefsPenalized,Ys,Zs,Xs,Xs_Te,Ys_Te,Zs_Te)
  gc()
  }
}
mean(AUCs)
mean(Concs)
Ys_Te <- split(Ytest, f = Ztest$Stage)
Resp <- do.call('rbind',Ys_Te)
Conc <- Est.Cval(cbind(Resp[,1],Resp[,2],LP_Tot), tau = tau,  nofit=T)$Dhat
AUC <- survivalROC.C(Stime = Resp[,1], status = Resp[,2], marker = LP_Tot, predict.time = RelTime)$AUC
citation("survivalROC")

summary(Ztest$Stage)
Ys_Te[[4]]
class(Zs[[1]])
