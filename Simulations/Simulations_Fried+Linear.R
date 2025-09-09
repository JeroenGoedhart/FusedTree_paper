setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
gc()
source('SideFunctions_Update.R')
source('RidgeFunctions_Update.R')
#libraries
library(VGAM)
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
library(mvtnorm)
library(Matrix)

##### simulations #####
Nsim <- 500
p <- 500
p_Clin <- 5
N <- 100 # or set to N = 100 for other setting
Ntest <- 5000

load("correlationmatrix.Rdata")
set.seed(1)
ids <- sample(1:ncol(CorrX), p, replace = F)
CorrX <- CorrX[ids,ids]
remove(ids)

set.seed(344)
Xtest <- mvtnorm::rmvnorm(Ntest, sigma = CorrX)
colnames(Xtest) <- paste0("x",seq(1,p))

g <- function(z) {
  15 * sin(pi * z[,1] * z[,2]) + 10 * (z[,3] - 0.5)^2 +
    2 * exp(z[,4]) + 2 * z[,5]
}

set.seed(344)
Ztest <- matrix(runif(Ntest*p_Clin,0,1),nrow = Ntest, ncol = p_Clin)
colnames(Ztest) <- paste0("z",seq(1,p_Clin))
Ztest <- data.frame(Ztest)
set.seed(4)
betas <- matrix(rlaplace(p, 0, 75/p), ncol = 1)
var(g(Ztest))
var(Xtest %*% betas)

signal <- var(g(Ztest) + Xtest %*% betas)

noise <- 5 ^ 2
R2 <- signal / (signal + noise)
set.seed(2)
Ytest <- g(Ztest) + Xtest %*% betas + rnorm(Ntest, 0, 5); var(Ytest)
Ytest <- Ytest[,1]
var(Ytest)

MSE_ridge <- c()
MSE_FusReg <- c()
MSE_FullFus <- c()
MSE_RF <- c()
MSE_GB <- c()
MSE_Lasso <- c()
MSE_ZeroFus <- c()

FusPar <- matrix(NA, nrow = Nsim, ncol = 5)
colnames(FusPar) <- c("Lambda","alpha","Lam Full Fus", "Lam Ridge","Lam Zero Fus")

for (i in 1:Nsim) {
  print(paste("simulation",i, sep = " "))
 
   set.seed(i ^ 3 + 534 + i)

  Z <- matrix(runif(N*p_Clin,0,1),nrow = N, ncol = p_Clin)
  colnames(Z) <- paste0("z",seq(1,p_Clin))
  Z <- data.frame(Z)
  X <- mvtnorm::rmvnorm(N, sigma = CorrX)
  colnames(X) <- paste0("x",seq(1,p))
  Y <- g(Z) + X %*% betas + rnorm(N,0,1)
  Y <- Y[,1]

  #specify folds for hyperparameter tuning
  folds <- CVfolds(Y = Y, model = "linear", balance = F, kfold = 10,nrepeat = 1)

  ####### 1. Fit Random Forest ######
  DFtrain <- data.frame(Ydf = Y,Xdf = cbind(X, Z))
  DFtest <- data.frame(Ydf = Ytest,Xdf = cbind(Xtest, Ztest))

  RF <- rfsrc(Ydf ~ .,data=DFtrain,ntree = 1000, mtry = floor((p+p_Clin)/3),
              var.used = "all.trees",importance = c("none"))
  preds_RF <- predict.rfsrc(RF, newdata = DFtest)$predicted
  MSE_RF[i] <- mean((preds_RF-Ytest)^2)
  remove(DFtrain,DFtest,RF,preds_RF)

  ####### 2. Fit gradient boosting ######
  tune_grid <- expand.grid(
    nrounds = c(100, 500, 1000),
    eta = c(0.01, 0.1, 0.3),
    max_depth = c(4),
    gamma = c(0),
    colsample_bytree = 0.75,
    min_child_weight = 30,
    subsample = 1
  )

  tune_control <- caret::trainControl(
    method = "cv", # cross-validation
    number = 3, # with n folds
    verboseIter = FALSE, # no training log
    allowParallel = FALSE #
  )

  xgb_tune <- caret::train(
    x = cbind(X,Z),
    y = Y,
    trControl = tune_control,
    tuneGrid = tune_grid,
    method = "xgbTree",verbosity = 0,verbose = F
  )
  xgb_tune$bestTune
  gbpred <- predict (xgb_tune,cbind(Xtest, Ztest))
  MSE_GB[i] <- mean((gbpred-Ytest)^2)
  remove(xgb_tune,gbpred,tune_control,tune_grid)


  ####### 3. Fit Ridge regression (Clinical unpenalized) ######
  Lam1 <- porridge::optPenaltyGLM.kCVauto(Y = Y, X = X, U = as.matrix(Z),
                                       lambdaInit = 10, lambdaGinit = 0, Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
                                       model = "linear", folds = folds, loss = "sos")

  if(Lam1 == Inf){
    Lam1 <- 10e8
  }

  RidgeFit <- ridgeGLM(Y = Y, U = as.matrix(Z), X = X,
                       lambda = Lam1, lambdaG = 0, Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
                       model = "linear")
  preds_Ridge <- (cbind(as.matrix(Ztest), Xtest) %*% RidgeFit)[,1]

  MSE_ridge[i] <- mean((preds_Ridge-Ytest)^2);
  remove(RidgeFit,preds_Ridge)

  ####### 2. Fit Lasso ######
  Las <- cv.glmnet(cbind(as.matrix(Z), X),Y, alpha = 1, nfolds = 10,family = "gaussian",
                  penalty.factor = c(rep(0, p_Clin), rep(1, p)))$lambda.min


  LassoFit <- glmnet(cbind(as.matrix(Z), X),Y, alpha = 1, lambda = Las, family = "gaussian",
                    penalty.factor = c(rep(0 ,p_Clin), rep(1, p)))

  preds_Lasso <- predict.glmnet(LassoFit, newx = cbind(as.matrix(Ztest), Xtest), s = Las, type = "response")[,1]

  MSE_Lasso[i] <- mean((preds_Lasso - Ytest)^2)
  remove(LassoFit,preds_Lasso)

  ####### 1. Fit Reg Tree ######
  dat <- cbind.data.frame(Y, Z)

  # fitting tree
  rp <- rpart(Y ~ ., data = dat, control = rpart.control(xval = 5, minbucket = 20),
              model = T)
  cp <- rp$cptable[,1][which.min(rp$cptable[,4])]
  Treefit <- prune(rp, cp = cp)
  remove(dat,rp)

  ####### 4. Fit fully FusedReg Tree ######
  Nodes <- row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
  Nodes <- as.numeric(Nodes)
  NumNod <- length(unique(Nodes))
  if (NumNod < 2){
    MSE_FullFus[i] <- MSE_ridge[i]
    MSE_ZeroFus[i] <- MSE_ridge[i]
    MSE_FusReg[i] <- MSE_ridge[i]
    FusPar[i,] <- c(c(Lam1, Inf), Lam1, Lam1, Lam1)
  } else {
    names <- paste0("N",sort(unique(Nodes)))

    Intercepts <- model.matrix(~0+factor(Nodes))
    colnames(Intercepts)[1:length(unique(Nodes))] <- names
    Intercepts <- as.matrix(cbind(Intercepts,Z))
    Lam <- porridge::optPenaltyGLM.kCVauto(Y = Y, X = X, U = Intercepts,
                                           lambdaInit = 100, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)),
                                           model = "linear", folds = folds, loss = "sos",
                                           maxIter = 100)



    Fit_Full <- ridgeGLM(Y = Y, U = Intercepts, X = X,
                         lambda = Lam, lambdaG = 0, Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
                         model = "linear")

    Nodes_test <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit, Ztest)]
    Nodes_test <- as.numeric(Nodes_test)
    Intercepts_Test <- model.matrix(~ 0 + factor(Nodes_test))
    Intercepts_Test <- as.matrix(cbind(Intercepts_Test, Ztest))
    colnames(Intercepts_Test)[1:length(unique(Nodes))] <- names


    Ypred_Full <- cbind(Intercepts_Test,Xtest) %*% Fit_Full
    MSE_FullFus[i] <- mean((Ypred_Full[,1]-Ytest)^2)

    remove(Intercepts,Intercepts_Test,Ypred_Full,Nodes,Nodes_test,Fit_Full)

    ####### 5. Fit  Zero FusedReg Tree ######
    Dat <- Dat_Tree(tree = Treefit, X = X, Z = Z, Y = Y, model = "linear",LinVars = T)
    X1 <- Dat$Omics; U1 <- Dat$Clinical; Y1 <- Dat$Response
    remove(Dat)
    optPenal <- porridge::optPenaltyGLM.kCVauto(Y = Y1, X = as.matrix(X1), U = U1,
                                                lambdaInit = 100, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X1), nrow = ncol(X1)),
                                                model="linear",
                                                folds=folds,
                                                loss="sos", maxIter=100)
    Fit_Zero <- porridge::ridgeGLM(Y = Y1, U = U1, X = as.matrix(X1),
                                   lambda = optPenal, lambdaG = 0, Dg = matrix(0, ncol=ncol(X1), nrow = ncol(X1)),
                                   model="linear")
    Dat_test <- Dat_Tree(tree = Treefit, X = Xtest, Z = Ztest, Y = Ytest,model = "linear",LinVars = T)
    X1_te <- Dat_test$Omics; U1_te <- Dat_test$Clinical; Y1_te = Dat_test$Response
    remove(Dat_test)

    Ypred_Zero <- cbind(U1_te,X1_te) %*% Fit_Zero
    MSE_ZeroFus[i] <- mean((Ypred_Zero[,1]-Ytest)^2)
    remove(U1,Y1,X1,X1_te,U1_te,Y1_te,Fit_Zero,Ypred_Zero)
    
    ####### 6. Fit  FusedTree Tree ######
    optPenalties <- PenOpt(Tree = Treefit, X = X, Y = Y, Z = Z, LinVars = T,
                           lambdaInit = Lam, alphaInit = 1e3,
                           model = "linear", folds = folds, loss = "sos")
    Fit <- FusTreeFit(Tree = Treefit, X = X, Y = Y, Z = Z, LinVars = T,
                     lambda = optPenalties[1], alpha = optPenalties[2],
                     model = "linear")


    Preds <- Predictions(fit = Fit, newX = Xtest, newZ = Ztest, newY = Ytest, Dat = T)
    MSE_FusReg[i] <- mean((Preds$Preds$Ypred-Preds$Preds$Resp)^2)
    remove(Preds,Fit)
    FusPar[i,] <- c(optPenalties, Lam,Lam1,optPenal)
    remove(optPenalties,optPenal,Lam,Lam1)

  }
  remove(X, Y, Z, Treefit, folds)
  gc()
}

results <- cbind.data.frame(MSE_FullFus,MSE_FusReg,MSE_GB,MSE_RF,MSE_ridge,MSE_Lasso,MSE_ZeroFus,VarY = rep(var(Ytest),Nsim), FusPar)
colMeans(results[1:7])
nm <- paste(N, Nsim, "FriedmanLinear_HighNoise.Rdata", sep = "_")
save(results, file = nm)

###### plotting results #####
library(ggplot2); library(viridis)

#1. boxplots of PMSE
load("100_500_FriedmanLinear_HighNoise.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1) <- c("FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
colMeans(results1)
dat <- stack(results1)
dat <- cbind.data.frame(dat, "Setting" = "N = 100")
load("300_500_FriedmanLinear_HighNoise.Rdata")
results1 <- results[,c(2,1,7,3,4,5,6)]
colnames(results1) <- c("FusTree","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
colMeans(results1)
dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1, "Setting" = "N = 300")

dat <- rbind.data.frame(dat, dat1)

name <- paste("PerformancePlot","FriedLin.pdf", sep = "_")
pdf(name, width = 7, height = 3)
bp <- ggplot(dat, aes(x = ind, y = values, group = ind)) +
  geom_boxplot(aes(fill = ind), outlier.shape = NA, fatten = 2) + coord_cartesian(ylim =  c(12, 45)) +
  scale_y_continuous(breaks = seq(12, 45, 8)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "D")+
  theme(legend.title = element_blank(), axis.text = element_text(size = 6)) + labs(x = "",y = "PMSE")
bp + facet_grid(. ~ Setting)
dev.off()

#2. diff PMSE versus alpha

load("100_500_FriedmanLinear_HighNoise.Rdata")
results$alpha[is.infinite(results$alpha)] <- 10e20
x <- log(results$alpha)
y <- results$MSE_FullFus/results$MSE_FusReg
dat <- cbind.data.frame("alpha" = x, "PMSE" = y, "Setting" = rep("N = 100", 500))


load("300_500_FriedmanLinear_HighNoise.Rdata")
max(results$alpha) #
x <- log(results$alpha)
y <- results$MSE_FullFus/results$MSE_FusReg
dat1 <- cbind.data.frame("alpha" = x, "PMSE" = y, "Setting" = rep("N = 300", 500))
dat <- rbind.data.frame(dat, dat1)

name <- paste("Perf_VS_alpha","FriedLin.pdf", sep = "_")
CairoPDF(name, width = 7, height = 3)
bp <- ggplot(dat, aes(x = alpha, y = PMSE)) +
  geom_point() +
  theme_light() +
  coord_cartesian(xlim =  c(-5, 15)) +
  theme(legend.title = element_blank()) + labs(x = "log(\u03b1)",y = "PMSE Ratio")
bp + facet_grid(. ~ Setting)
dev.off()
