setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/FusedTree_Paper/Applications")

#BiocManager::install("mcsurvdata")
library(mcsurvdata)
library(ExperimentHub)
library(mice) #some missing values
library(survival)
library(viridis)

####################################
######## Data Preprocessing ########
####################################

# load data in environment using mcsurvdata package
eh <- ExperimentHub()
nda.crc <- query(eh, "mcsurvdata")[["EH1498"]]
dat<-nda.crc
remove(eh,nda.crc)
gc()

### define Omics
Omics <- t(exprs(dat))
dim(Omics)
Features <- featureData(dat)
colnames(Omics)<- Features@data[["Gene.Symbol"]]
length(unique(colnames(Omics)))
Omics <- Omics[,-which(duplicated(colnames(Omics)))]
dim(Omics)
gc()

### define Clinical
Clinical <- pData(phenoData(dat))
Clinical<-cbind.data.frame(gender = Clinical$gender,
                           age = Clinical$age,
                           site = Clinical$site,
                           Stage = Clinical$stage,
                           CMS = Clinical$cms,
                           time = Clinical$tev,
                           event = Clinical$evn)

# remove missing response values
ids = which(is.na(Clinical$time) | is.na(Clinical$event))
Clinical <- Clinical[-ids,]; Omics <- Omics[-ids,]
dim(Omics)

## impute missing values with single imputation
Response = Surv(time = Clinical$time, event = Clinical$event)
Response[,1][which(Response[,1]==0)] <- 0.0001
Clinical <- Clinical[,-c(6,7)]

set.seed(1)
Imp = mice(data = Clinical, m = 1)
Clinical<-complete(Imp)
noNArows<-apply(is.na(Clinical),1,sum)==0
noNAcols<-apply(is.na(Clinical),2,sum)==0
Clinical$gender<-as.numeric(Clinical$gender)- 1 #0 = male; 1 = female
Clinical$site <- as.numeric(Clinical$site) - 1# 0 = right; 1 = left-rectum
Clinical$Stage = factor(Clinical$Stage, levels(Clinical$Stage))

save(Clinical,Omics,Response, file = "CRCdata.Rdata")


######################################
######## Preliminary Analysis ########
######################################
setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/FusedTree_Paper/Applications")
load("CRCdata.Rdata")
library(survival)
library(survminer)

ClinFit <- coxph(Response ~.,data = Clinical)
summary(ClinFit)


fit <- survfit(Response ~ 1, data = Clinical)
sum(fit$n.event)
name <- paste("SurvCurve_CRC.pdf",sep = "_")
pdf(name, width=3.25,height=3.2)
plot1 <- ggsurvplot(fit, data = Clinical, palette = "#440154FF")$plot
plot1
dev.off()

#######################################
##############  Analysis ##############
#######################################
setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/FusedTree_Paper/Applications")
load("CRCdata.Rdata")

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


# split data in training and test set
set.seed(48)
ids=sample(1:nrow(Clinical),size = 0.2*nrow(Clinical),replace = F)
Y<- Surv(time = Response[,1], event = Response[,2])
X<-Omics[-ids,]; Xtest<-Omics[ids,]
Y<-Response[-ids,]; Ytest<-Response[ids,]
Z<-Clinical[-ids,]; Ztest=Clinical[ids,]
remove(ids,Response,Clinical,Omics)
gc()
dim(X)

# Fit training and test survival curve
d.demo = data.frame(1)
fit <- survfit(Y ~ 1, data = d.demo)
fit1 <- survfit(Ytest ~ 1, data = d.demo)
sum(fit$n.event)
sum(fit1$n.event)
fitTot <- list(Training = fit, Test = fit1)

name <- paste("SurvCurve_CRC.pdf",sep = "_")
pdf(name, width=3.25,height=3, onefile = F)
ggsurvplot(fitTot, data = d.demo, combine = TRUE, # Combine curves
           risk.table = F,                  # Add risk table
           conf.int = TRUE,                    # Add confidence interval
           conf.int.style = "step",            # CI style, use "step" or "ribbon"
           censor = FALSE,                     # Remove censor points
           tables.theme = theme_cleantable(),  # Clean risk table
           palette = "uchicago",
           legend.labs=c("Train", "Test"),
           legend.title="",
           legend = c(0.2,0.2),
           font.x = c(10, face = "bold"),
           font.y = c(10, face = "bold"))
dev.off()



##### Model Fitting #####
#########################

source('Sidefunctions_Update.R')
source('RidgeFunctions_Update.R')


p_clin = ncol(model.matrix(~.,Z)[,-1])
p = ncol(X)


Lin = T
RelTime = 5
tau = 8


# 0. Clinical cox ph
ClinFit <- coxph(Y ~.,data = Z)
summary(ClinFit)
LP <- predict(ClinFit,newdata = Z,type = "lp")
LpPred <- predict(ClinFit,newdata = Ztest,type = "lp")
#ConcClin <- concordance(Ytest ~ LpPred)$concordance
ConcClin <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LpPred), tau = tau,  nofit=T)$Dhat
AUCClin <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LpPred, predict.time = RelTime)$AUC
remove(ClinFit,LpPred,LP)

#1. RF Clin + Omics
DFtrain <- cbind.data.frame(time=Y[,1], event=Y[,2],X,Z)
DFtest <- cbind.data.frame(time=Ytest[,1],event=Ytest[,2],Xtest,Ztest)
set.seed(4)
p <- ncol(X)
CV <- tune(Surv(time,event) ~ .,data=DFtrain, ntreeTry = 500, mtryStart = p/20,
           nodesizeTry = c(5,10,20,100))
nsize <- CV$optimal[1]
mtry <- CV$optimal[2]
set.seed(4)
RF <- rfsrc(Surv(time,event) ~ .,data=DFtrain, ntree = 500, var.used="all.trees",importance=c("none"),splitrule="logrank",
            mtry = mtry, nodesize = nsize)
preds_RF <- predict.rfsrc(RF, newdata = DFtest, outcome = "train")
ConcRF <-  Est.Cval(cbind(Ytest[,1],Ytest[,2],rowSums(preds_RF$chf)), tau = tau, nofit=T)$Dhat
AUCRF <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = rowSums(preds_RF$chf), predict.time = RelTime)$AUC
remove(DFtrain,DFtest,preds_RF,RF,CV, mtry,nsize)


#2. Ridge
set.seed(1)
foldsHyp <- CVfolds(Y=Y, model = "cox", kfold = 5, nrepeat = 3)
Lam1 <- optPenaltyGLM.kCVauto2(Y = Y, X = X, U = model.matrix(~0+., Z), 
                            lambdaInit = 100, 
                            lambdaGinit = 0, maxIter = 100, minSuccDiff = 1e-5,
                            model = "surv", folds = foldsHyp, loss = "loglik")
RidgeFit <- ridgeGLM2(Y = Y, U = model.matrix(~ 0 + ., Z), X = X, 
                      lambda = Lam1, lambdaG = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)), 
                      model="surv")
names(RidgeFit) <- c(colnames(model.matrix(~0+., Z)),colnames(X))
LP = (cbind(model.matrix(~0+., Ztest),Xtest) %*% RidgeFit)[,1]
ConcRidge <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LP), tau = tau, nofit=T)$Dhat
AUCRidge <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP, predict.time = RelTime)$AUC
remove(RidgeFit,LP)
gc()
#3. Residual Approach 
ClinFit <- coxph(Y ~.,data = Z)
summary(ClinFit)
LP <- predict(ClinFit,newdata = Z,type = "lp")
LpPred <- predict(ClinFit,newdata = Ztest,type = "lp")
Las1 <- cv.glmnet(X,Y, offset = matrix(LP, ncol = 1, nrow = nrow(X)), alpha = 0, nfolds = 5,family = "cox")$lambda.min
ResFit <- glmnet(X,Y, offset = matrix(LP, ncol = 1, nrow = nrow(X)),  alpha = 0, lambda = Las1, family = "cox")
LP1 <- LpPred + (Xtest %*% ResFit$beta[,1])[,1]
Conc_Res <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LP1), tau = tau, nofit=T)$Dhat
AUC_Res <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP1, predict.time = RelTime)$AUC
remove(ResFit,LP,LP1)

#4. Lasso
Las = cv.glmnet(cbind(model.matrix(~.,Z)[,-1],X),Y, alpha = 1, nfolds = 5,family = "cox",
                standardize = TRUE,penalty.factor = c(rep(0,p_clin),rep(1,p)))$lambda.min
LassoFit = glmnet(cbind(model.matrix(~.,Z)[,-1],X),Y, alpha = 1, lambda = Las, family = "cox",
                  standardize = TRUE, penalty.factor = c(rep(0,p_clin),rep(1,p)))

LP = (cbind(model.matrix(~., Ztest)[,-1],Xtest) %*% LassoFit$beta)[,1]
Conc_Lasso <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LP), tau = tau, nofit=T)$Dhat
AUC_Lasso <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP, predict.time = RelTime)$AUC
remove(LassoFit,LP)
gc()

#5. boosting #
max.depths <- c(2,4,6)
etas <- c(0.1, 0.01, 0.001)
Res <- expand.grid(Depth = max.depths, Shrinkage = etas)
Res <- cbind(Res, CV.Perf = NA) # add column for cv performances
NTrees <- 100

for (i in 1:nrow(Res)) {
  print(paste("fold ", i))
  depth <- Res$Depth[i]
  eta <- Res$Shrinkage[i]
  set.seed(i)
  perf <- c()
  for (j in 1:length(foldsHyp)){
    gbm <- gbm::gbm.fit(x = cbind.data.frame(Z,X), y = Y,
                        n.trees = NTrees,
                        interaction.depth = depth,
                        shrinkage = eta,
                        distribution="coxph",
                        nTrain = 541, verbose = F)
    perf[j] <- gbm$valid.error[NTrees]
  }
  Res[i,3] <- mean(perf[j])
}
id <- which.max(Res$CV.Perf)
#id <- which.min(Res$CV.Perf)
Res[id,]

set.seed(11)
gbm <- gbm.fit(x = cbind.data.frame(Z,X), y = Y,
               n.trees = NTrees,
               interaction.depth = Res$Depth[id],
               shrinkage = Res$Shrinkage[id],
               distribution="coxph", verbose = F)
gbm <- gbm.fit(x = cbind.data.frame(Z,X), y = Y,
               n.trees = NTrees,
               interaction.depth = 2,
               shrinkage = 0.001,
               distribution="coxph", verbose = F)
set.seed(11)
gbm.pred <- predict.gbm(gbm, 
                        newdata=cbind.data.frame(Ztest,Xtest), 
                        type="response") 

ConcGB <- Est.Cval(cbind(Ytest[,1],Ytest[,2],gbm.pred), tau = tau, nofit=T)$Dhat
AUCGB <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = gbm.pred, predict.time = RelTime)$AUC
remove(gbm,gbm.pred,Res)

#### BlockForest ####
#####################
gc()

Xblock = cbind.data.frame(Z,X)
XblockTest = cbind.data.frame(Ztest,Xtest)
colnames(Xblock) <- paste("X", 1:ncol(Xblock), sep="")
colnames(XblockTest) <- paste("X", 1:ncol(XblockTest), sep="")

blocks <- rep(1:2, times=c(ncol(Z), ncol(X)))
blocks <- lapply(1:2, function(x) which(blocks==x))
set.seed(40)
blockforobj <- blockfor(Xblock, Y,  blocks=blocks,
                        nsets = 50)
blockforobj$paramvalues
set.seed(40)
Preds=predict(blockforobj$forest, data=XblockTest)
Preds <- rowSums(Preds$chf)
ConcBlock <- Est.Cval(cbind(Ytest[,1],Ytest[,2],Preds), tau = tau, nofit=T)$Dhat
AUCBlock <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = Preds, predict.time = RelTime)$AUC
remove(Xblock,XblockTest,blocks,,Preds,blockforobj)
gc()

Inf.Cval.Delta(cbind(Ytest[,1],Ytest[,2]), Preds, LP, tau = tau, itr = 10000, seed = NULL)

##### FusedTree ####
####################

dat=cbind.data.frame(Y,Z)
set.seed(9)
rp <- rpart(Y~.,data = dat, control = rpart.control(xval =5, minbucket = 30,cp=0),
            model = T)
minerr <- which.min(rp$cptable[,"xerror"])
bestcp <- rp$cptable[minerr,"CP"]
remove(minerr)
Treefit <- prune(rp, cp = bestcp)
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
foldsHyp = .CVfoldsTree(Y=Y,Tree = Treefit,Z=Z,model="surv", kfold = 5, nrepeat = 3)
optLam1 = PenOpt(Tree=Treefit,Y=Y,X=X,Z=Z,model = "surv",lambdaInit = Lam1/NumNod, alphaInit = 100, 
                 folds = foldsHyp,LinVars = Lin, maxIter = 30)
optLam1
Fit <- FusTreeFit(Tree=Treefit,Y=Y,X=as.matrix(X),Z=Z,
                  model = "surv",
                  lambda = optLam1[1],
                  alpha =  optLam1[2],
                  LinVars = Lin)

Preds=Predictions(fit=Fit,newX = as.matrix(Xtest),newZ = Ztest,newY = Ytest, model="surv", Linvars = Lin)
Preds = Preds$Preds
lpmin=as.numeric(Preds$LinPred)

ConcFusTree=Est.Cval(cbind(Ytest[,1],Ytest[,2],lpmin), tau=tau, nofit=T)$Dhat
AUCFusTree <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = lpmin, predict.time = RelTime)$AUC

### Fully Fused ###
Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
Nodes <- as.numeric(Nodes)
NumNod <- length(unique(Nodes))
Intercepts = model.matrix(~0+factor(Nodes))
Intercepts <- as.matrix(cbind(Intercepts,Z$age))
LamFully <-optPenaltyGLM.kCVauto2(Y = Y, X = X, U = Intercepts, 
                                       lambdaInit = 100, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)),
                                       model ="surv", folds = foldsHyp, loss = "loglik")
LamFully


Fit_Full <- ridgeGLM2(Y = Y, U = Intercepts, X = X, 
                     lambda = LamFully, lambdaG = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)), 
                     model="surv")

Nodes_test <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Ztest)]
Nodes_test <- as.numeric(Nodes_test)
Intercepts_Test = model.matrix(~0+factor(Nodes_test))
Intercepts_Test <- as.matrix(cbind(Intercepts_Test,Ztest$age))
Ypred_Full = cbind(Intercepts_Test,Xtest) %*% Fit_Full
ConcFully <- Est.Cval(cbind(Ytest[,1],Ytest[,2],Ypred_Full), tau=tau, nofit=T)$Dhat
AUCFully <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = Ypred_Full, predict.time = RelTime)$AUC
remove(Intercepts,Intercepts_Test,Ypred_Full,Nodes,Nodes_test,Fit_Full,Ypred_Full)


### backward ###
Nodenames = paste0("N",sort(unique(Nodes)))
pvals = sapply(sort(unique(Nodes)), function(x) globaltest::p.value(globaltest::gt(
  Surv(Y[Nodes==x,1],Y[Nodes==x,2]), 
  alternative = data.frame(matrix(X[Nodes == x, ],nrow =length(which(Nodes==x)),ncol=ncol(X))),
  model = "cox")))
names(pvals) <- Nodenames
pvals
pvals[is.na(pvals)]<-1
print("Estimated p values in the nodes are:")
print(pvals)
print("Remove omics effects in nodes in order:")
pvals<-sort(pvals,decreasing = T)
print(paste(names(pvals),collapse = " -> "))
Dat = Dat_Tree(tree = Treefit, X = X, Y = Y,Z = Z,  LinVars = Lin,model = "surv")
X1 = Dat$Omics; Y1 = Dat$Response; U1 = Dat$Clinical
colSums(U1)
Dat_tree_test = Dat_Tree(tree=Treefit, X=Xtest, Y=Ytest, Z=Ztest, LinVars = Lin,model = "surv")
Ute_tree = Dat_tree_test$Clinical; Yte_tree=Dat_tree_test$Response; Xte_tree = Dat_tree_test$Omics

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
    if (ncol(Delta) != ncol(X2)){stop("number of columns of penalty matrix does not equal number of columns of design matrix X")}
    #if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
    #if(alphaInit <= 0){stop("initial penalty alphaInit should be larger than zero")}
    optPenalties <- optPenaltyGLM.kCVauto2(Y = Y1, X = X2, U = U1, 
                                           lambdaInit = 1500, lambdaGinit = 10000,
                                           folds = foldsHyp,
                                           Dg=Delta,model = "surv",loss="loglik")
    
    print(optPenalties)
    Fit = ridgeGLM2(Y = Y1, U = U1, X = X2, 
                    lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                    model="surv")
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
optPenalties
nm <- "FinalFitCRC.Rdata"
save(Fit,X2,U1,Y, file = nm)

# Fully Fused #



###########################################################
############### Downstream Analysis #######################
###########################################################

load("FinalFitCRC.Rdata")
betas = Fit[-(1:7)]
betas
beta_N5 <- betas[seq(1,length(betas),3)]
beta_N12 <- betas[seq(2,length(betas),3)]
beta_N13 <- betas[seq(3,length(betas),3)]
sum(abs(beta_N5))
sum(abs(beta_N12))
sum(abs(beta_N13))

## gene signature
signature <- c("ESCO2_N", "AXIN2_N", "PLK1_N", "CDC25C_N", "IGF1_N", "TREX2_N", "ALKBH2_N", "ESR1_N",  "MC1R_N")
ids <- c()

for (i in 1:length(signature)){
  nm <- signature[i]
  id = which(grepl(nm,names(Fit)))
  print(names(Fit)[id])
  ids <- append(ids,id)
}
names(Fit)[ids]

pathway <- Fit[ids]
pathway
N5 <- pathway[seq(1,length(pathway),3)]
N12 <- pathway[seq(2,length(pathway),3)]   
N13 <- pathway[seq(3,length(pathway),3)]
abs(sum(N5))
abs(sum(N12))
sum(abs(N13))

############
beta_Node <- split(betas, ceiling(seq_along(betas)/3))
beta_Node[[1]]
unlist()
vars <- unlist(lapply(beta_Node, sd))

ord <- order(vars,decreasing = T)
vars[tail(ord)]
beta_Node[ord[20:40]]
beta_Node[ord[21292]]

##### transcrips
Trans = c("ZEB","SNAI1","SNAI2","TWIST1")
ids <- c()

for (i in 1:length(Trans)){
  nm <- Trans[i]
  id = which(grepl(nm,names(Fit)))
  print(names(Fit)[id])
  ids <- append(ids,id)
}
pathway <- Fit[ids]
pathway
N5 <- pathway[seq(1,length(pathway),3)]
N12 <- pathway[seq(2,length(pathway),3)]   
N13 <- pathway[seq(3,length(pathway),3)]
abs(sum(N5))
abs(sum(N12))
sum(abs(N13))

##### reg paths #####
library(reshape2)
genes <- c("MAGEA6","HLA-DRB4")
id1 = which(grepl("MAGEA6",names(RidgeFit)))
names(RidgeFit)
ids <- c()
ids1 <- c()
for (i in 1:length(genes)){
  nm <- genes[i]
  id = which(grepl(nm,names(Fit)))
  #id1 = which(grepl(nm,names(RidgeFit)))
  print(names(Fit)[id])
  ids <- append(ids,id)
  #ids1 <- append(ids1,id1)
}
RidgeFit[ids1]
Lam = 1508
alf = alpha = seq(0.001,10e5,500)
Betas = lapply(alf, function (x) ridgeGLM2(Y = Y1, U = U1, X = X2, 
                                           lambda = Lam, lambdaG = x, Dg = Delta, 
                                           model="surv"))
Betas = do.call(rbind,Betas)
colnames(Betas)= c(colnames(U1),colnames(X2))
save(Betas, file = "Betas_RegPath_CRC.Rdata")
load("Betas_RegPath_CRC.Rdata")
Betas <- Betas[,ids]

colnames(Betas) <- paste0("\u03b2_", colnames(Betas))
colnames(Betas) <- gsub('N','',colnames(Betas))

Res = cbind.data.frame("Alpha"=log(alf),Betas)
melt = melt(Res, id = "Alpha")


name <- paste("RegPath","CRC", sep = "_")
library(ggplot2); library(viridis)
library(Cairo)
CairoPDF(name, width=3.25,height=2.9)
plot <- ggplot(melt, aes(x = Alpha, y = value, group = variable, colour = variable,linetype = variable)) + 
  geom_line(linewidth=1.1) + theme_light() +
  scale_color_viridis(discrete = T, option = "D") + scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3))) +
  theme(legend.title = element_blank(),legend.text=element_text(size=5)) + labs(x="log(\u03b1)",y="Effect estimates") 
plot <- plot + geom_vline(xintercept=log(14836), linetype="dotted") 
plot
dev.off()
library(grid)
vp.BottomRight <- viewport(height=unit(.9, "npc"), width=unit(0.5, "npc"), 
                           just=c("left","top"),
                           y=.9, x=0.5)
library(Cairo)
# plot your base graphics 

name <- paste("Plot_Application_CRC.pdf")
CairoPDF(name, width=8,height=3)
par(mfrow=c(1,2))
rpart.plot(Treefit, # middle graph
           type=5,
           extra=1, 
           box.palette="Pu",
           branch.lty=8, 
           shadow.col=0, 
           nn=TRUE,
           cex = 0.6)


# plot the ggplot using the print command
print(plot, vp=vp.BottomRight)
mtext(side = 2, line = 2, "a", cex = 1, font = 2,las = 2 ,at = 1)
mtext(side = 4, line = 2, "b", cex = 1, font = 2,las = 2 ,at = 1)
mtext("X", side = 1, line = 3.7, col = "red", at=0.05, cex = 1.2)
mtext("X", side = 1, line = 3.7, col = "red", at=0.22, cex = 1.2)
mtext("X", side = 1, line = 3.7, col = "red", at=0.941, cex = 1.2)
dev.off()
getwd()

##### Stratified ####
#####################

summary(Z$Stage)
Stage = levels(Z$Stage)
Resp <- c()
LP <- c()

for (i in 1:length(Stage)) {
  subset <- Stage[i]
  ids = which(Z$Stage==subset)
  ids1 = which(Ztest$Stage==subset)
  X1 <- X[ids,]; Y1 <- Y[ids,]
  X1_Te <- Xtest[ids1,]; Y1_Te <- Ytest[ids1,]
  
  if (i==1){
    Lam <- 1
  } else {
    Lam <- penalized::optL2(response = Y1, penalized = X1,
                            model = "cox", fold = 10)
    Lam <- Lam$lambda
  }
  RidgeFit1 = penalized::penalized(response = Y1, penalized = X1,
                                   model = "cox", lambda2=1)
  
  coefsPenalized <- penalized::coefficients(RidgeFit1)
  LpPred <- (X1_Te %*% coefsPenalized)[,1]
  LP <- append(LP,LpPred)
  Resp <- append(Resp,Y1_Te)
}
LP <- -LP
ConcStrat <- Est.Cval(cbind(Resp[,1],Resp[,2],LP), tau = tau,  nofit=T)$Dhat
AUCStrat <- survivalROC.C(Stime = Resp[,1], status = Resp[,2], marker = LP, predict.time = RelTime)$AUC


