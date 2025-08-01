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
library(glinternet)

##### simulations #####
Nsim <- 500
p <- 500
p_Clin <- 5
N <- 300 # or set to N = 100 for other setting
Ntest <- 5000

load("correlationmatrix.Rdata")
set.seed(1)
ids <- sample(1:ncol(CorrX), p, replace = F)
CorrX <- CorrX[ids,ids]
remove(ids)

set.seed(344)
Xtest <- mvtnorm::rmvnorm(Ntest, sigma = CorrX)
colnames(Xtest) <- paste0("x",seq(1, p))
set.seed(344)
Ztest <- matrix(runif(Ntest * p_Clin, 0, 1),nrow = Ntest, ncol = p_Clin)
Ztest[,1] <- sample(0:3, Ntest, replace = T)
Ztest[,2] <- sample(0:1, Ntest, replace = T)
colnames(Ztest) <- paste0("z",seq(1, p_Clin))
Ztest <- data.frame(Ztest)
set.seed(4)
betas <- matrix(rlaplace(p, 0, 30/p),ncol = 1)

g <- function (z,x,betas){
  1*(z[,1] < 1.5)*(z[,2] < 0.5)*(-10 + x[,1:(p/4)] %*% betas[1:(p/4)]*8) + #node1
    1*(z[,1] < 1.5)*(z[,2] >= 0.5)*(-5 + x[,1:(p/4)]%*% betas[1:(p/4)]*2) + #node2
    1*(z[,1] >= 1.5)*(z[,3] < 0.5)*(5 + x[,1:(p/4)]%*% betas[1:(p/4)]*0.5) +
    1*(z[,1] >= 1.5)*(z[,3] >= 0.5)*(10 + x[,1:(p/4)]%*% betas[1:(p/4)]*0.125) +
    x[,((p/4)+1):p]%*% betas[((p/4)+1):p] + 5*z[,4]
}
set.seed(4)
Ytest <- g(z = Ztest, x = Xtest, betas = betas)[,1] + rnorm(Ntest,0,1)
var(Ytest)
full_signal <- g(Ztest, Xtest, betas)
var_signal <- var(full_signal)
true_R2 <- var_signal / (var_signal + 1)
# Z-only intercept component (as previously computed)
intercept_only <- 
  1 * (Ztest[,1] < 1.5) * (Ztest[,2] < 0.5) * (-10) +
  1 * (Ztest[,1] < 1.5) * (Ztest[,2] >= 0.5) * (-5) +
  1 * (Ztest[,1] >= 1.5) * (Ztest[,3] < 0.5) * (5) +
  1 * (Ztest[,1] >= 1.5) * (Ztest[,3] >= 0.5) * (10)

# Additive Z part
Z_additive <- 5 * Ztest[,4]

# Pure X contribution = signal - Z-only parts
X_contribution <- full_signal - intercept_only - Z_additive

# Variance explained by X
var_X_contribution <- var(X_contribution)
var(Z_additive+intercept_only)


treef <- function(z){
  nodes <- ifelse((z[,1] < 1.5 & z[,2] < 0.5),1,0) +
    ifelse((z[,1] < 1.5 & z[,2] >= 0.5),2,0) +
    ifelse((z[,1] >= 1.5 & z[,3] < 0.5),3,0) +
    ifelse((z[,1] >= 1.5 & z[,3] >= 0.5),4,0)
  return(nodes)
}

NodesTe_oracle <- treef(Ztest)
NodesTe_oracle <- factor(NodesTe_oracle, levels = sort(unique(NodesTe_oracle)))
summary(NodesTe_oracle)
Ute_or <- model.matrix(~ 0 + NodesTe_oracle)
Ute_or <- cbind(Ute_or, Ztest[,-c(1:3)])
Ute_or <- as.matrix(Ute_or)
remove(NodesTe_oracle)

MSE_ridge <- c()
MSE_FusReg <- c()
MSE_FullFus <- c()
MSE_RF <- c()
MSE_GB <- c()
MSE_Lasso <- c()
MSE_ZeroFus <- c()
MSE_oracle <- c()
MSE_glinternet <- c()
AllNodes <- c()
FusPar <- matrix(NA,nrow = Nsim, ncol = 7)
colnames(FusPar) <- c("Lambda","alpha","Lam Full Fus", "Lam Ridge","Lam Zero Fus","Lambda_Or","alpha_Or")
VarsBeta <- matrix(NA,nrow = Nsim, ncol = 2)

for (i in 1:Nsim) {
  
  print(paste("simulation",i, sep = " "))
  set.seed(i^3 + 5343)
 
  Z <- matrix(runif(N * p_Clin, 0, 1),nrow = N, ncol = p_Clin) 
  Z[,1] <- sample(0:3, N, replace = T)
  Z[,2] <- sample(0:1, N, replace = T)
  colnames(Z) <- paste0("z", seq(1, p_Clin))
  Z <- data.frame(Z)
  X <- mvtnorm::rmvnorm(N, sigma = CorrX)
  colnames(X) <- paste0("x", seq(1, p))
  
  ## simulating response
  Y <- g(z = Z, x = X,betas = betas) + rnorm(N, 0, 1)
  Y <- Y[,1]
  
  #specify folds for hyperparameter tuning
  folds <- CVfolds(Y = Y, model = "linear", balance = F, kfold = 10,nrepeat = 1)
  
  ####### 1. Fit Random Forest ######
  DFtrain <- data.frame(Ydf = Y, Xdf = cbind(X, Z))
  DFtest <- data.frame(Ydf = Ytest, Xdf = cbind(Xtest, Ztest))
  
  RF <- rfsrc(Ydf ~ .,data = DFtrain, ntree = 500, mtry = floor((p + p_Clin) / 3),
              var.used = "all.trees", importance = c("none"))
  preds_RF <- predict.rfsrc(RF, newdata = DFtest)$predicted
  MSE_RF[i] <- mean((preds_RF - Ytest)^2)
  remove(DFtrain, DFtest, RF, preds_RF)
  
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
    x = cbind(X, Z),
    y = Y,
    trControl = tune_control,
    tuneGrid = tune_grid,
    method = "xgbTree",verbosity = 0,verbose = FALSE
  )
  gbpred <- predict (xgb_tune, cbind(Xtest, Ztest))
  MSE_GB[i] <- mean((gbpred - Ytest)^2)
  remove(xgb_tune, gbpred, tune_control, tune_grid)
  
  
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
  preds_Ridge <- (cbind(as.matrix(Ztest),Xtest) %*% RidgeFit)[,1]
  
  MSE_ridge[i] <- mean((preds_Ridge - Ytest)^2)
  remove(RidgeFit, preds_Ridge)
  
  ####### 4. Fit Lasso ######
  Las <- cv.glmnet(cbind(as.matrix(Z), X), Y, alpha = 1, nfolds = 10,family = "gaussian",
                  penalty.factor = c(rep(0, p_Clin), rep(1, p)))$lambda.min
  
  
  LassoFit <- glmnet(cbind(as.matrix(Z), X), Y, alpha = 1, lambda = Las, family = "gaussian",
                    penalty.factor = c(rep(0, p_Clin), rep(1, p)))
  
  preds_Lasso <- predict.glmnet(LassoFit, newx = cbind(as.matrix(Ztest),Xtest), s = Las, type = "response")[,1]
  
  MSE_Lasso[i] <- mean((preds_Lasso - Ytest)^2)
  remove(LassoFit, preds_Lasso)
  
  ####### 4. Fit glinternet ######
  numLevs <- c(4, 2, 1, 1, 1, rep(1, ncol(X)))
  Lam <- glinternet.cv(X = cbind(as.matrix(Z), X), Y = Y, numLevels = numLevs, nFolds = 5,
                       interactionCandidates = c(1:5), family = "gaussian")
  
  fit <- glinternet(X = cbind(as.matrix(Z), X), Y = Y, numLevels = numLevs, lambda = Lam$lambdaHat,
                    interactionCandidates = c(1:5), family = "gaussian")
  Ypred <- predict(fit, X = cbind(as.matrix(Ztest), Xtest), lambda = Lam$lambdaHat)
  MSE_glinternet[i] <- mean((Ypred - Ytest)^2)
  remove(Lam, fit, Ypred)
  
  #### oracle Tree Model ####
  Nodes_oracle <- treef(Z)
  Nodes_oracle <- factor(Nodes_oracle, levels = sort(unique(Nodes_oracle)))
  summary(Nodes_oracle)
  U_or <- model.matrix(~ 0 + Nodes_oracle)
  X_or <- t(Matrix::KhatriRao(t(X), t(U_or)))
  U_or <- cbind(U_or, Z[,-c(1, 2, 3)])
  U_or <- as.matrix(U_or)
  Delta_or <- .PenMatr(NumNodes = length(unique(Nodes_oracle)), p = p)
  
  
  optPenalties_or <- optPenaltyGLM.kCVauto2(Y = Y, X = as.matrix(X_or), U = U_or, 
                                            lambdaInit = 10, lambdaGinit = 10, Dg=as.matrix(Delta_or),
                                            model = "linear",
                                            folds = folds,
                                            loss = "sos")
  
  
  Fit_or <- ridgeGLM2(Y = Y, U = U_or, X = as.matrix(X_or), 
                      lambda = optPenalties_or[1], lambdaG = optPenalties_or[2], Dg = as.matrix(Delta_or), 
                      model = "linear")
  
  Xte_Or <- t(Matrix::KhatriRao(t(Xtest), t(Ute_or[,1:4])))
  Ypred <- cbind(Ute_or, as.matrix(Xte_Or)) %*% as.matrix(Fit_or)
  Ypred <- Ypred[,1]
  MSE_oracle[i] <- mean((Ypred - Ytest)^2)
  remove(Xte_Or,X_or,Delta_or,Ypred,Nodes_oracle)
  
  ####### 1. Fit Reg Tree ######
  dat <- cbind.data.frame(Y, Z)
  # fitting tree
  rp <- rpart(Y ~ .,data = dat, control = rpart.control(xval = 5,  maxdepth = 2),
              model = T)
  cp <- rp$cptable[,1][which.min(rp$cptable[,4])]
  Treefit <- prune(rp, cp = cp)
  rpart.plot(Treefit, # middle graph
             type = 5,
             extra = 1, 
             box.palette = "Pu",
             branch.lty = 8, 
             shadow.col = "gray", 
             nn = TRUE,
             cex = 0.6)
  remove(dat,rp)
  
  ####### 4. Fit fully FusedReg Tree ######
  Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit, Z)]
  Nodes <- as.numeric(Nodes)
  NumNod <- length(unique(Nodes))
  AllNodes[i] <- NumNod
  if (NumNod < 2){
    MSE_FullFus[i] <- MSE_ridge[i]
    MSE_ZeroFus[i] <- MSE_ridge[i]
    MSE_FusReg[i] <- MSE_ridge[i]
    FusPar[i,] <- c(c(Lam1, Inf), Lam1, Lam1, Lam1, optPenalties_or)
  } else {
    names <- paste0("N",sort(unique(Nodes)))
    Intercepts <- model.matrix(~ 0 + factor(Nodes))
    colnames(Intercepts)[1:length(unique(Nodes))] <- names
    Intercepts <- as.matrix(cbind(Intercepts, Z[,-1]))
    Lam <- porridge::optPenaltyGLM.kCVauto(Y = Y, X = X, U = Intercepts, 
                                           lambdaInit = 10, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)),
                                           model = "linear", folds = folds, loss = "sos")
    Fit_Full <- ridgeGLM(Y = Y, U = Intercepts, X = X, 
                         lambda = 10, lambdaG = 0, Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)), 
                         model = "linear")
    
    Nodes_test <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Ztest)]
    Nodes_test <- as.numeric(Nodes_test)
    Intercepts_Test <- model.matrix(~ 0 + factor(Nodes_test))
    Intercepts_Test <- as.matrix(cbind(Intercepts_Test, Ztest[,-1]))
    colnames(Intercepts_Test)[1:length(unique(Nodes))] <- names
    
    Ypred_Full <- cbind(Intercepts_Test,Xtest) %*% Fit_Full
    MSE_FullFus[i] <- mean((Ypred_Full[,1] - Ytest)^2)
    remove(Intercepts,Intercepts_Test,Ypred_Full,Nodes,Nodes_test,Fit_Full)
    
    ####### 5. Fit  Zero FusedReg Tree ######
    Dat <- Dat_Tree(tree = Treefit, X = X, Z = Z, Y = Y, model = "linear",LinVars = T)
    X1 <- Dat$Omics; U1 <- Dat$Clinical; Y1 <- Dat$Response
    remove(Dat)
    optPenal <- porridge::optPenaltyGLM.kCVauto(Y = Y1, X = as.matrix(X1), U = U1, 
                                                lambdaInit = 10, lambdaGinit = 0, Dg = matrix(0, ncol = ncol(X1), nrow = ncol(X1)),
                                                model ="linear",
                                                folds = folds,
                                                loss = "sos")
    
    Fit_Zero <- porridge::ridgeGLM(Y = Y1, U = U1, X = as.matrix(X1), 
                                   lambda = optPenal, lambdaG = 0, Dg = matrix(0, ncol = ncol(X1), nrow = ncol(X1)), 
                                   model = "linear")
    Dat_test <- Dat_Tree(tree = Treefit, X = Xtest, Z = Ztest, Y = Ytest,model = "linear", LinVars = T)
    X1_te <- Dat_test$Omics; U1_te <- Dat_test$Clinical; Y1_te <- Dat_test$Response
    remove(Dat_test)
    Ypred_Zero <- cbind(U1_te, X1_te) %*% Fit_Zero
    MSE_ZeroFus[i] = mean((Ypred_Zero[,1] - Ytest)^2)
    remove(U1,Y1,X1,X1_te,U1_te,Y1_te,Fit_Zero,Ypred_Zero)
    
    ####### 6. Fit  FusedTree ######
    optPenalties <- PenOpt(Tree = Treefit, X = X, Y = Y, Z = Z, LinVars = T,
                           lambdaInit = 100, alphaInit = 10,
                           model = "linear", folds = folds, loss = "sos")
    
    Fit <- FusTreeFit(Tree = Treefit, X = X, Y = Y, Z = Z, LinVars = T,
                     lambda = optPenalties[1], alpha = optPenalties[2],
                     model = "linear")
    
    Preds <- Predictions(fit = Fit, newX = Xtest, newZ = Ztest, newY = Ytest, Dat = T)
    
    Fit1 <- Fit$Effects[-seq(1,ncol(Z)+length(unique(Fit$Tree$where)))]
    a <- split(Fit1, ceiling(seq_along(Fit1)/NumNod))
    variances <- unlist(lapply(a,var))
    VarsBeta[i,] <- c(mean(variances[1:(p/4)]),mean(variances[((p/4)+1):p]))
    remove(Fit1,a,variances)
    MSE_FusReg[i] <- mean((Preds$Preds$Ypred - Preds$Preds$Resp)^2)
    remove(Preds,Fit)
    FusPar[i,] <- c(optPenalties,  Lam, Lam1, optPenal, optPenalties_or)
    remove(optPenalties,optPenal,Lam,Lam1)
  }
  
  remove(X,Y,Z,Treefit,folds)
  gc()
}

mean(MSE_FullFus); mean(MSE_FusReg); mean(MSE_glinternet); mean(MSE_ZeroFus); mean(MSE_GB); mean(MSE_oracle)
save(AllNodes,file = "Nodes_300.Rdata")

results <- cbind.data.frame(MSE_FullFus,MSE_FusReg,MSE_GB,MSE_RF,MSE_ridge,MSE_Lasso,MSE_oracle,MSE_ZeroFus, MSE_glinternet,FusPar,VarsBeta)
results <- cbind.data.frame(MSE_FullFus,MSE_FusReg,MSE_GB,MSE_RF,MSE_ridge,MSE_Lasso,MSE_oracle,MSE_ZeroFus, MSE_glinternet)
colMeans(results)[1:10]
nm <- paste(N, Nsim, "EffModd_New.Rdata", sep = "_")
save(results, file = nm)

###### plotting results #####
library(ggplot2); library(viridis)

#1. boxplots of PMSE
load("100_500_EffModd1.Rdata")
results1 <- results[,c(2,7,1,8,3,4,6)]
colnames(results1)<- c("FusTree","Oracle","FullFus","ZeroFus","GB","RF","Lasso")


dat <- stack(results1)
dat <- cbind.data.frame(dat, "Setting" = "N = 100")
load("300_500_EffModd1.Rdata")

results1 <- results[,c(2,7,1,8,3,4,6)]
colnames(results1)<- c("FusTree","Oracle","FullFus","ZeroFus","GB","RF","Lasso")
colMeans(results1)
mean(results$MSE_ridge)

dat1 <- stack(results1)
dat1 <- cbind.data.frame(dat1, "Setting" = "N = 300")

dat <- rbind.data.frame(dat, dat1)


name <- paste("PerformancePlot", "EffMod.pdf", sep = "_")
pdf(name, width = 7, height = 3)
bp <- ggplot(dat, aes(x = ind, y = values, group = ind)) +
  geom_boxplot(aes(fill = ind), outlier.shape = NA,fatten = 2) + coord_cartesian(ylim =  c(10, 42.5)) +
  scale_y_continuous(breaks = seq(10, 42.5, 5)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "D")+
  theme(legend.title = element_blank(), axis.text = element_text(size = 7)) + labs(x = "",y = "PMSE") 
bp + facet_grid(. ~ Setting)
dev.off()

#2. diff PMSE versus alpha
load("100_500_EffModd1.Rdata")
x <- log(results$alpha)

y <- results$MSE_ZeroFus/results$MSE_FusReg
dat <- cbind.data.frame("alpha" = x, "PMSE" = y,"Setting" = rep("N = 100",500))
load("300_500_EffModd1.Rdata")
x <- log(results$alpha)
y <- results$MSE_ZeroFus/results$MSE_FusReg
dat1 <- cbind.data.frame("alpha" = x, "PMSE" = y,"Setting" = rep("N = 300",500))
dat <- rbind.data.frame(dat, dat1)


name <- paste("Perf_VS_alpha","EffMod.pdf", sep = "_")
CairoPDF(name, width = 7, height = 3)
bp <- ggplot(dat, aes(x = alpha, y = PMSE)) +
  geom_point() + 
  theme_light() +
  coord_cartesian(xlim =  c(-10, 15)) +
  scale_x_continuous(breaks = seq(-10, 15, 5)) +
  theme(legend.title = element_blank()) + labs(x = "log(\u03b1)", y = "PMSE Ratio") 
bp + facet_grid(. ~ Setting)
dev.off()


load("100_500_EffModd1.Rdata")
colMeans(results)
load("300_500_EffModd.Rdata")
colMeans(results)
plot(log(results$Lambda),log(results$`Lam Zero Fus`))
x <- log(results$alpha)
y <- results$MSE_ZeroFus/results$MSE_FusReg
dat <- cbind.data.frame("lambda" = x, "PMSE" = y,"Setting" = rep("N = 100", 500))

name <- paste(300, 500, "Perf_VS_alpha","EffMod.pdf", sep = "_")
pdf(name, width = 7, height = 3)
bp <- ggplot(dat, aes(x = lambda, y = PMSE)) +
  geom_point() +
  theme_light() +
  theme(legend.title = element_blank()) + labs(x = "log(\u03b1)", y = "PMSE Ratio") 
bp
dev.off()

load("300_500_EffModd.Rdata")
results1 <- results[,c(2,7,1,8,3,4,5,6)]
colnames(results1)<- c("FusTree","Oracle","FullFus","ZeroFus","GB","RF","Ridge","Lasso")
dat <- stack(results1)
name <- paste(300, 500, "BoxplotPerf","EffMod.pdf", sep = "_")
pdf(name, width = 7, height = 3)
bp <- ggplot(dat, aes(x = ind, y = values)) +
  geom_boxplot(aes(fill = ind), outlier.shape = NA, fatten = 2) + coord_cartesian(ylim =  c(10, 35)) +
  scale_y_continuous(breaks=seq(10,35,5)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "D", begin = 0.1) +
  theme(legend.title = element_blank(), axis.text = element_text(size = 6)) + labs(x = "", y = "PMSE") 
bp
dev.off()

load("100_500_EffModd.Rdata")
res <- cbind.data.frame("\u03bb FusedTree" = log(results$Lambda),
                        "\u03bb ZeroFus" = log(results$`Lam Zero Fus`))
res<- stack(res)
res <- cbind.data.frame(res,"Setting" = "N = 100")
load("300_500_EffModd.Rdata")
res1 <- cbind.data.frame("\u03bb FusedTree" = log(results$Lambda),
                        "\u03bb ZeroFus" = log(results$`Lam Zero Fus`))
res1<- stack(res1)
res1 <- cbind.data.frame(res1,"Setting" = "N = 300")

results = rbind.data.frame(res, res1)


bp <- ggplot(results, aes(x = values, fill = ind)) +
  geom_histogram(bins = 50) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "D", begin = 0.5, end = 0.9) +
  theme(legend.title = element_blank()) + labs(x = "log (\u03bb)", y = "Frequency") 
bp + facet_grid(. ~ Setting)