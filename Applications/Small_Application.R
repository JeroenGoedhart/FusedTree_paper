### Application 2 ###

setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Applications")
#install.packages("EBcoBART")
library(EBcoBART)
library(glmnet)
library(randomForestSRC)
library(porridge)
library(rpart); 
library(rpart.plot)
library(treeClust)
library(multiridge)
library(ModTools)
library(caret)
library(xgboost)
library(glinternet)

library(fusedTree)
library(globaltest)

data(Lymphoma)
Ytest <- Lymphoma$Ytest; Y <- Lymphoma$Ytrain
Xtrain <- Lymphoma$Xtrain; Xtest <- Lymphoma$Xtest
remove(Lymphoma) ; gc()

Z <- Xtrain[,137:140]; Ztest <- Xtest[,137:140]
X <- as.matrix(Xtrain[,-c(137:140)]); Xtest <- as.matrix(Xtest[,-c(137:140)])
remove(Xtrain)
#Z$IPI <- factor(Z$IPI, ordered = T)
#Ztest$IPI <- factor(Ztest$IPI, ordered = T)
#Ztest$IPI[Ztest$IPI == 0] <- 1 # as zero does not happen for training data
### 1. Random Forest (classification) ###
DFtrain <- data.frame(Ydf = factor(Y), Xdf = cbind(X, Z))
DFtest  <- data.frame(Ydf = factor(Ytest), Xdf = cbind(Xtest, Ztest))


p <- ncol(X)
q <- ncol(Z)
set.seed(1)
CV <- tune(Ydf ~ .,data=DFtrain, ntreeTry = 1000, mtryStart = p/10,
           nodesizeTry = c(5,10,20,100))
nsize <- CV$optimal[1]
mtry <- CV$optimal[2]
RF <- rfsrc(Ydf ~ ., data = DFtrain, ntree = 1000, mtry = mtry,
            nodesize = nsize,
            var.used = "all.trees", importance = "none", family = "class")

preds_RF_prob <- predict.rfsrc(RF, newdata = DFtest)$predicted[, 2]  # probability for class 1
AUC_RF <- auc(Ytest, preds_RF_prob)
Brier_RF <- mean((preds_RF_prob - Ytest)^2)
remove(RF, preds_RF_prob, DFtrain, DFtest,mtry,nsize,CV)

### 3. Ridge Regression (binomial, clinical unpenalized) ###
set.seed(2)
folds <- CVfolds(Y=Y, model = "logistic", balance = F, kfold = 5, nrepeat = 3)
Lam1 <- porridge::optPenaltyGLM.kCVauto(
  Y = Y, X = X, U = model.matrix(~.,Z),
  lambdaInit = 10, lambdaGinit = 0,
  Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
  model = "logistic", folds = folds)

if (Lam1 == Inf) Lam1 <- 1e8

RidgeFit <- ridgeGLM(
  Y = Y, U = model.matrix(~.,Z), X = X,
  lambda = Lam1, lambdaG = 0,
  Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
  model = "logistic"
)

ridge_prob <- (cbind(model.matrix(~.,Ztest), Xtest) %*% RidgeFit)[, 1]
ridge_prob <- 1 / (1 + exp(-ridge_prob))  # convert logit to probability
AUC_ridge <- auc(Ytest, ridge_prob)
Brier_ridge <- mean((ridge_prob - Ytest)^2)
remove(RidgeFit, ridge_prob,folds)

### 4. Lasso (logistic) ###
Las <- cv.glmnet(cbind(model.matrix(~.,Z), X), Y, alpha = 1, nfolds = 5, family = "binomial",
                 penalty.factor = c(rep(0, 5), rep(1, p)))$lambda.min

LassoFit <- glmnet(cbind(model.matrix(~.,Z), X), Y, alpha = 1, lambda = Las, family = "binomial",
                   penalty.factor = c(rep(0, 5), rep(1, p)))

lasso_prob <- predict(LassoFit, newx = cbind(model.matrix(~.,Ztest), as.matrix(Xtest)),
                      s = Las, type = "response")[, 1]
AUC_Lasso <- auc(Ytest, lasso_prob)
Brier_Lasso <- mean((lasso_prob - Ytest)^2)
remove(LassoFit, lasso_prob)


####### 6. Fit fusedTree ######

dat=cbind.data.frame(Y,Z)
testdata=cbind.data.frame(Ytest,Ztest)
dat$Y <- as.factor(dat$Y)
testdata$Y <- as.factor(testdata$Y)
# fitting tree
Treefit <- rpart(Y~.,data = dat)
name <- paste("TreeLymphoma.pdf",sep = "_")
pdf(name, width=5,height=3, onefile = F)
rpart.plot(Treefit, # middle graph
           type=5,
           extra=1, 
           box.palette="Pu",
           branch.lty=8, 
           shadow.col=0, 
           nn=TRUE,
           cex = 0.6)
dev.off()

Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
Nodes <- as.numeric(Nodes)
NumNod <- length(unique(Nodes))
Nodenames = paste0("N",sort(unique(Nodes)))
pvals = sapply(sort(unique(Nodes)), function(x) globaltest::p.value(globaltest::gt(
  Y[Nodes==x], 
  alternative = data.frame(matrix(X[Nodes == x, ],nrow =length(which(Nodes==x)),ncol=ncol(X))),
  model = "logistic")))
names(pvals) <- Nodenames
pvals



probs <- predict(Treefit, newdata = testdata, type = "prob")[, 2]
AUC_tree <- as.numeric(pROC::auc(Ytest, probs))
Brier_tree <- mean((probs-Ytest)^2)
remove(dat)
set.seed(30)
foldsTree <- CVfoldsTree(Y = Y, Tree = Treefit, Z = data.frame(Z), model = "logistic", nrepeat = 3, kfold = 5)

# 6a. Regular fusedTree
optPenalties <- PenOpt(Treefit, X, Y, data.frame(Z), "logistic", lambdaInit = 1000, alphaInit = 10,
                       loss = "loglik", LinVars = FALSE, folds = foldsTree, multistart = T)
optPenalties
fit <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "logistic",
                 lambda = optPenalties[1], alpha = optPenalties[2], verbose = FALSE)

Ypred_prob <- predict(fit, newX = Xtest, newY = Ytest, newZ = Ztest)$Ypred
AUC_fusedTree <- as.numeric(pROC::auc(Ytest, Ypred_prob))
Brier_fusedTree <- mean((Ytest - Ypred_prob)^2)

#### Clinical only ####
### 1. Random Forest (classification) ###
DFtrain <- data.frame(Ydf = factor(Y), Xdf = Z)
DFtest  <- data.frame(Ydf = factor(Ytest), Xdf = Ztest)
set.seed(8)
CV <- tune(Ydf ~ .,data=DFtrain, ntreeTry = 500, mtryStart = 3,
           nodesizeTry = c(3,5,10,20,100))
nsize <- CV$optimal[1]
mtry <- CV$optimal[2]
RF <- rfsrc(Ydf ~ ., data = DFtrain, ntree = 500, mtry = mtry,
            nodesize = nsize,
            var.used = "all.trees", importance = "none", family = "class")

preds_RF_prob <- predict.rfsrc(RF, newdata = DFtest)$predicted[, 2]  # probability for class 1
AUC_RF1 <- auc(Ytest, preds_RF_prob)
Brier_RF1 <- mean((preds_RF_prob - Ytest)^2)
remove(RF, preds_RF_prob, DFtrain, DFtest,mtry,nsize)

Xblock = cbind.data.frame(Z,X)
XblockTest = cbind.data.frame(Ztest,Xtest)
colnames(Xblock) <- paste("X", 1:ncol(Xblock), sep="")
colnames(XblockTest) <- paste("X", 1:ncol(XblockTest), sep="")

blocks <- rep(1:2, times=c(ncol(Z), ncol(X)))
blocks <- lapply(1:2, function(x) which(blocks==x))
library(blockForest)
set.seed(40)
blockforobj <- blockfor(Xblock, Y,  blocks=blocks)
blockforobj1 <- blockfor(Xblock, Y,  blocks=blocks,
                         nsets = 100, always.select.block = 1)
blockforobj$paramvalues
set.seed(40)
Preds=predict(blockforobj1$forest, data=XblockTest)

Preds <- Preds$predictions
Brier_Block <- mean((Preds-Ytest)^2)
AUC_Block <- as.numeric(auc(Ytest,Preds))
remove(Xblock,XblockTest,blocks,Preds,blockforobj)




### 2. Gradient Boosting (XGBoost) ###
### 2. Gradient Boosting (XGBoost) ###
Ztest$IPI <- as.numeric(Ztest$IPI)
Ztest <- as.matrix(Ztest)
tune_grid <- expand.grid(
  nrounds = c(100, 200,500),
  eta = c(0.01, 0.05, 0.1, 0.3),
  max_depth = c(2,3,4),
  gamma = 0,
  colsample_bytree = 0.75,
  min_child_weight = c(2,5,8,10),
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = FALSE
)

Y_class <- factor(ifelse(Y == 1, "yes", "no"))  # required format for caret
xgb_tune <- caret::train(
  x = cbind(Z,X),
  y = Y_class,
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  verbosity = 0,
  metric = "ROC"
)
xgb_tune$bestTune

xgb_pred_prob <- predict(xgb_tune, cbind(Ztest, Xtest), type = "prob")[, "yes"]
AUC_GB <- auc(Ytest, xgb_pred_prob)
Brier_GB <- mean((xgb_pred_prob - Ytest)^2)



tune_grid <- expand.grid(
  nrounds = c(100, 200,500),
  eta = c(0.01,0.05, 0.1, 0.3),
  max_depth = c(2,3),
  gamma = 0,
  colsample_bytree = 0.75,
  min_child_weight = c(5,10,15),
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = FALSE
)

Y_class <- factor(ifelse(Y == 1, "yes", "no"))  # required format for caret
xgb_tune <- caret::train(
  x = Z,
  y = Y_class,
  trControl = tune_control,
  tuneGrid = tune_grid,
  method = "xgbTree",
  verbosity = 0,
  metric = "ROC"
)
xgb_tune$bestTune

xgb_pred_prob <- predict(xgb_tune, Ztest, type = "prob")[, "yes"]
AUC_GB1 <- auc(Ytest, xgb_pred_prob)
Brier_GB1 <- mean((xgb_pred_prob - Ytest)^2)
remove(xgb_tune, xgb_pred_prob)

### glm ###
dat=cbind.data.frame(Y,Z)
testdata=cbind.data.frame(Ytest,Ztest)
dat$Y <- as.factor(dat$Y)
testdata$Y <- as.factor(testdata$Y)

# Fit logistic regression model
glm_fit <- glm(Y ~ ., data = dat, family = binomial)

# Predict probabilities
probs <- predict(glm_fit, newdata = testdata, type = "response")

# AUC
library(pROC)

AUC_glm <- auc(Ytest, probs)
print(auc$auc)

# Brier Score
Brier_glm <- mean((probs - Ytest)^2)
print(brier_score)


