
setwd("C:/Users/VNOB-0732/Desktop/R files/ClinOmics/Applications")
gc()
#BiocManager::install("mcsurvdata")
library(mcsurvdata)
library(ExperimentHub)
library(mice) #some missing values
library(survival)
library(viridis)

####################################
######## Data Preprocessing ########
####################################

eh <- ExperimentHub()
nda.crc <- query(eh, "mcsurvdata")[["EH1498"]]
#nda.brca <- query(eh, "mcsurvdata")[["EH1497"]]
dat<-nda.crc
#dat<-nda.brca
remove(eh,nda.crc)
#remove(eh,nda.brca)
gc()

### Omics
Omics <- t(exprs(dat))
dim(Omics)
Features <- featureData(dat)
colnames(Omics)<- Features@data[["Gene.Symbol"]]
length(unique(colnames(Omics)))
Omics<- Omics[,-which(duplicated(colnames(Omics)))]
dim(Omics)

############
### Clinical
Clinical <- pData(phenoData(dat))
summary(Clinical$cms)
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
sum(noNArows)
sum(noNAcols)
Clinical$gender<-as.numeric(Clinical$gender)- 1 #0 = male; 1 = female
Clinical$site <- as.numeric(Clinical$site) - 1# 0 = right; 1 = left-rectum
Clinical$Stage = factor(Clinical$Stage, levels(Clinical$Stage))

Omics1 <- Omics
#p_sel <- 500
#p_sel <- 5000
p_sel <- ncol(Omics1) #only use p_sel omics covariates with largest variance
sds = apply(Omics1,2,function(x) sd(x,na.rm = T))
ord = order(sds,decreasing = T)[1:p_sel]
Omics <- Omics1[,ord]
remove(sds,ord)

nm <- paste(p_sel,"CRC.Rdata", sep = "_")
dim(Omics)
save(Clinical,Omics,Response, file = nm)
remove(nm)
remove(p_sel)

######################################
######## Preliminary Analysis ########
######################################
load("500_CRC.Rdata")
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
load("500_CRC.Rdata")
load("5000_CRC.Rdata")
load("21292_CRC.Rdata")
dim(Omics)
gc()

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

# standardize data
Omics <- scale(Omics)
#set.seed(4)

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
# 0. Clinical
ClinFit <- coxph(Y ~.,data = Z)
summary(ClinFit)
LP = predict(ClinFit,newdata = Z,type = "lp")
LpPred = predict(ClinFit,newdata = Ztest,type = "lp")
#ConcClin <- concordance(Ytest ~ LpPred)$concordance
ConcClin1 <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LpPred), tau = tau,  nofit=T)$Dhat
AUCClin <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LpPred, predict.time = RelTime)$AUC

remove(ClinFit,LpPred)
#1. RF
DFtrain <- cbind.data.frame(time=Y[,1], event=Y[,2],X,Z)
DFtest <- cbind.data.frame(time=Ytest[,1],event=Ytest[,2],Xtest,Ztest)
RF <- rfsrc(Surv(time,event) ~ .,data=DFtrain,ntree=2000, var.used="all.trees",importance=c("none"),splitrule="logrank",
            mtry = 4*sqrt(p))
preds_RF <- predict.rfsrc(RF, newdata = DFtest, outcome = "train")

#ConcRF <- 1-preds_RF$err.rate[2000]
ConcRF1 <-  Est.Cval(cbind(Ytest[,1],Ytest[,2],rowSums(preds_RF$chf)), tau = tau, nofit=T)$Dhat
AUCRF <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = rowSums(preds_RF$chf), predict.time = RelTime)$AUC
remove(DFtrain,DFtest,preds_RF,RF)

#1. RF
DFtrain <- cbind.data.frame(time=Y[,1], event=Y[,2],X)
DFtest <- cbind.data.frame(time=Ytest[,1],event=Ytest[,2],Xtest)
RF <- rfsrc(Surv(time,event) ~ .,data=DFtrain,ntree=2000, var.used="all.trees",importance=c("none"),splitrule="logrank",
            mtry = 4*sqrt(p))
preds_RF <- predict.rfsrc(RF, newdata = DFtest, outcome = "train")

ConcRF <- 1-preds_RF$err.rate[2000]
ConcRF1 <-  Est.Cval(cbind(Ytest[,1],Ytest[,2],rowSums(preds_RF$chf)), tau = tau, nofit=T)$Dhat
AUCRF <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = rowSums(preds_RF$chf), predict.time = RelTime)$AUC
remove(DFtrain,DFtest,preds_RF,RF)


#2. Ridge
set.seed(1)
foldsHyp <- CVfolds(Y=Y, model = "cox", kfold = 5,nrepeat = 3)
start=proc.time()
Lam1=optPenaltyGLM.kCVauto2(Y = Y, X = X, U = model.matrix(~0+., Z), 
                            lambdaInit = 100, 
                            lambdaGinit = 0, maxIter = 100, minSuccDiff = 1e-5,
                            model = "surv", folds = foldsHyp, loss = "loglik")

end=proc.time()-start
end
#start=proc.time()
#Lam = penalized::optL2(response = Y, penalized = X, unpenalized = model.matrix(~0+., Z)[,-1],
#            model = "cox", lambda1 = 100, fold = 5)
#end1=proc.time()-start
#end
#end1

#Lam <- Lam$lambda
#RidgeFit1 = penalized::penalized(response = Y, penalized = X,
#                                unpenalized = model.matrix(~., Z)[,-1],
#                                model = "cox", lambda2=Lam1)

#coefsPenalized <- penalized::coefficients(RidgeFit1)


RidgeFit <- ridgeGLM2(Y = Y, U = model.matrix(~0+., Z), X = X, 
                      lambda = Lam*3, lambdaG = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)), 
                      model="surv")
gc()
Lam1
names(RidgeFit) <- c(colnames(model.matrix(~0+., Z)),colnames(X))
LP = (cbind(model.matrix(~0+., Ztest),Xtest) %*% RidgeFit)[,1]
LP1 = -LP

#ConcRidge <- concordance(Ytest ~ LP1)$concordance
ConcRidge1 <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LP), tau = tau, nofit=T)$Dhat
AUCRidge <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP, predict.time = RelTime)$AUC


Las = cv.glmnet(cbind(model.matrix(~.,Z)[,-1],X),Y, alpha = 1, nfolds = 5,family = "cox",
                penalty.factor = c(rep(0,p_clin),rep(1,p)))$lambda.min
Las1 = cv.glmnet(X,Y, offset = matrix(LP, ncol = 1, nrow = nrow(X)), alpha = 0, nfolds = 5,family = "cox")$lambda.min


LassoFit = glmnet(cbind(model.matrix(~.,Z)[,-1],X),Y, alpha = 1, lambda = Las, family = "cox",
                  penalty.factor = c(rep(0,p_clin),rep(1,p)))
LassoFit1 = glmnet(X,Y, offset = matrix(LP, ncol = 1, nrow = nrow(X)),  alpha = 0, lambda = Las1, family = "cox")
length(LassoFit1$beta)
LP = (cbind(model.matrix(~., Ztest)[,-1],Xtest) %*% LassoFit1$beta)[,1]
LP1 = LpPred + (Xtest %*% LassoFit1$beta[,1])[,1]
#ConcLasso <- concordance(Ytest ~ LP)$concordance
ConcLasso1 <- Est.Cval(cbind(Ytest[,1],Ytest[,2],LP1), tau = tau, nofit=T)$Dhat
AUCLasso <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = LP1, predict.time = RelTime)$AUC
remove(LassoFit,LP)
gc()
### boosting ###
DFtrain <- cbind.data.frame(time=Y[,1], event=Y[,2],X,Z)
DFtest <- cbind.data.frame(time=Ytest[,1],event=Ytest[,2],Xtest,Ztest)
gc()

set.seed(11)
options(expressions = 5e5)
gbm <- gbm(Surv(time, event)~.,
           DFtrain,
           n.trees = 100,
           interaction.depth=2,
           distribution="coxph")

set.seed(11)
gbm.pred <- predict.gbm(gbm, 
                        newdata=DFtest, 
                        type="response") 
gbm.pred <- -gbm.pred


#ConcGB <- concordance(Ytest ~ gbm.pred)$concordance
ConcGB1 <- Est.Cval(cbind(Ytest[,1],Ytest[,2],-gbm.pred), tau = tau, nofit=T)$Dhat
AUCGB <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = -gbm.pred, predict.time = RelTime)$AUC
remove(DFtrain,DFtest,gbm,gbm.pred)

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
Preds <- -rowSums(Preds$chf)
#ConcBlock <- concordance(Ytest ~ Preds)$concordance
ConcBlock1 <- Est.Cval(cbind(Ytest[,1],Ytest[,2],-Preds), tau = tau, nofit=T)$Dhat
AUCBlock <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = -Preds, predict.time = RelTime)$AUC
remove(Xblock,XblockTest,blocks,Preds,blockforobj)
gc()

##### CoxBoost ####
####################




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

name <- paste("Tree_CRC.pdf",sep = "_")
pdf(name, width=3.25,height=3)
rpart.plot(Treefit, # middle graph
           type=5,
           extra=1, 
           box.palette="Pu",
           branch.lty=8, 
           shadow.col="gray", 
           nn=TRUE,
           cex = 0.6)


dev.off()

Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
Nodes <- as.numeric(Nodes)
NumNod <- length(unique(Nodes))

#Lam1=1000
set.seed(5)
foldsHyp = .CVfoldsTree(Y=Y,Tree = Treefit,Z=Z,model="surv",kfold = 5, nrepeat = 3)
start = proc.time()
optLam1 = PenOpt(Tree=Treefit,Y=Y,X=X,Z=Z,model = "surv",lambdaInit = Lam1/NumNod,alphaInit = 1000, 
                 folds = foldsHyp,LinVars = Lin, maxIter = 30)
end1 = proc.time()-start
end1
optLam1
Fit=FusTreeFit(Tree=Treefit,Y=Y,X=as.matrix(X),Z=Z,model = "surv",lambda = optLam1[1],alpha =  optLam1[2],LinVars = Lin)

#View(as.matrix(Fit$Omics))
#View(Fit$Clinical)

Preds=Predictions(fit=Fit,newX = as.matrix(Xtest),newZ = Ztest,newY = Ytest, model="surv", Linvars = Lin)
Preds = Preds$Preds
lpmin=as.numeric(Preds$LinPred)

#ConcFusTree=concordance(Ytest ~ lpmin)$concordance
ConcFusTree1=Est.Cval(cbind(Ytest[,1],Ytest[,2],lpmin), tau=tau, nofit=T)$Dhat
AUCFusTree <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = lpmin, predict.time = RelTime)$AUC
Fit$Effects[1:9]
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
  j=3
  EmpNodes = names(pvals)[1:j]
  #names(Fits)[i] <- paste(names(pvals)[1:i],collapse = ",")
  print(paste("Fit FusedTree without omics effects in", paste(names(pvals)[1:j],collapse = ", ")))
  ids = lapply(EmpNodes, function (x) grepl(x,str_sub(colnames(X1), start = -2)))
  ids1 = which(Reduce("|",ids))
  
  X2 <- X1[,-ids1]
  Xte <- Xte_tree[,-ids1]
  remove(ids,ids1)
  NodNum <- ncol(X2)/p
  dim(X2)
  gc()
  if (NodNum == 0){
    Fit <- glm(Y1 ~ 0 + ., data = data.frame(U1),family = "binomial")$coefficients
    
    ids = which(names(Fit) %in% colnames(Ute_tree))
    Ests = Fit[ids]
    LP <- as.numeric(Ute_tree %*% Ests)
    Preds <- exp(LP)/(1+exp(LP))
    AUCs[j] <- roc(Yte_tree,Preds)$auc
    MSEs[j] <- mean((Preds-Yte_tree)^2)
    remove(Fit,ids)
    
  }
  if (NodNum == 1){
    dim(X2)
    
    Lam <- optPenaltyGLM.kCVauto2(Y = Y1, X = X2, U = U1, 
                                  lambdaInit = Lam1, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X2), nrow = ncol(X2)),
                                  folds = foldsHyp,
                                  model = "surv", loss="loglik")
    Fit <-ridgeGLM2(Y = Y1, U = U1, X = X2, 
                    lambda = Lam1, lambdaG = 0, Dg = matrix(0, ncol=ncol(X2), nrow = ncol(X2)), 
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
                                           lambdaInit = Lam1/NodNum, lambdaGinit = 1000,
                                           folds = foldsHyp,
                                           Dg=Delta,model = "surv",loss="loglik")
    
    optPenalties
    optLam1
    Fit = ridgeGLM2(Y = Y1, U = U1, X = X2, 
                    lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                    model="surv")
    names(Fit)<-c(colnames(U1),colnames(X2))
    LP = -as.numeric(cbind(Ute_tree,Xte) %*% Fit)
    Concordances1[j] <- Est.Cval(cbind(Ytest[,1],Ytest[,2],-LP), tau=tau, nofit=T)$Dhat
    AUCs[j] <- survivalROC.C(Stime = Ytest[,1], status = Ytest[,2], marker = -LP, predict.time = RelTime)$AUC
    remove(Fit)
  }
  remove(X2)
}
Concordances1
AUCs
optPenalties
nm <- "FinalFitCRC.Rdata"
save(Fit,X2,U1,Y, file = nm)

###########################################################
############### Downstream Analysis #######################
###########################################################

load("FinalFitCRC.Rdata")
betas = Fit[-(1:7)]
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
  id1 = which(grepl(nm,names(RidgeFit)))
  print(names(Fit)[id])
  ids <- append(ids,id)
  ids1 <- append(ids1,id1)
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


substring(colnames(Betas),1,1)[-1] <-"\u03b2"
Res = cbind.data.frame("Alpha"=log(alf),Betas)
melt = melt(Res, id = "Alpha")


name <- paste("RegPath","CRC", sep = "_")
library(ggplot2); library(viridis)
library(Cairo)
CairoPDF(name, width=3.25,height=2.9)
plot <- ggplot(melt, aes(x = Alpha, y = value, group = variable, colour = variable,linetype = variable)) + 
  geom_line() + theme_light() +
  scale_color_viridis(discrete = T, option = "D") + scale_linetype_manual(values = c(rep("solid",3),rep("dashed",3))) +
  theme(legend.title = element_blank(),legend.text=element_text(size=4)) + labs(x="log(\u03b1)",y="Effect estimates") 
plot <- plot + geom_vline(xintercept=log(14836), linetype="dotted") 
plot
dev.off()
