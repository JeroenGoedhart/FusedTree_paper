setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
gc()

#libraries
library(VGAM)
library(glmnet)
library(randomForestSRC)
library(porridge)
#library(penalized)
#library(magic)
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
library(fusedTree)
##### simulations #####

Nsim = 500
p = 500
p_Clin = 5
N = 300
Ntest = 5000
gc()
load("correlationmatrix.Rdata")
set.seed(1)
ids=sample(1:ncol(CorrX),p,replace = F)
CorrX <- CorrX[ids,ids]
remove(ids)



set.seed(344)
Xtest = mvtnorm::rmvnorm(Ntest, sigma = CorrX)
colnames(Xtest) = paste0("x",seq(1,p))
set.seed(344)
Ztest = matrix(runif(Ntest*p_Clin,0,1),nrow = Ntest, ncol = p_Clin)
Ztest[,1] <- sample(0:3,Ntest,replace = T)
Ztest[,2] <- sample(0:1,Ntest,replace = T)
colnames(Ztest) = paste0("z",seq(1,p_Clin))
Ztest = data.frame(Ztest)
set.seed(4)
betas <- matrix(rlaplace(p,0,10/p),ncol = 1)
#Var(Xtest[,1:25]%*%betas[1:25])

g <- function (z,x,betas){
  1*(z[,1] < 1.5)*(z[,2] < 0.5)*(-2 + x[,1:p/4] %*% betas[1:p/4]*8) + #node1
    1*(z[,1] < 1.5)*(z[,2] >= 0.5)*(-1 + x[,1:p/4]%*% betas[1:p/4]*2) + #node2
    1*(z[,1] >= 1.5)*(z[,3] < 0.5)*(1 + x[,1:p/4]%*% betas[1:p/4]*0.5) +
    1*(z[,1] >= 1.5)*(z[,3] >= 0.5)*(2 + x[,1:p/4]%*% betas[1:p/4]*0.125) +
    x[,(p/4+1):p]%*% betas[(p/4+1):p] 
}

logit_prob <- function(z, x, betas) {
  eta <- g(z = z, x = x, betas = betas)
  prob <- 1 / (1 + exp(-eta))  # sigmoid function
  return(prob)
}
set.seed(4)
prob_test <- logit_prob(z = Ztest, x = Xtest, betas = betas)
Ytest <- rbinom(Ntest, size = 1, prob = prob_test)
summary(Ytest)

# storage containers #
AUC_RF <- numeric(Nsim)
Brier_RF <- numeric(Nsim)

AUC_GB <- numeric(Nsim)
Brier_GB <- numeric(Nsim)

AUC_ridge <- numeric(Nsim)
Brier_ridge <- numeric(Nsim)

AUC_Lasso <- numeric(Nsim)
Brier_Lasso <- numeric(Nsim)

AUC_glinternet <- numeric(Nsim)
Brier_glinternet <- numeric(Nsim)

AUC_fusedTree <- numeric(Nsim)
Brier_fusedTree <- numeric(Nsim)

AUC_fulFus <- numeric(Nsim)
Brier_fulFus <- numeric(Nsim)

AUC_zeroFus <- numeric(Nsim)
Brier_zeroFus <- numeric(Nsim)


n_leaves_vec <- numeric(Nsim)
alphas <- numeric(Nsim)

for (i in 1:Nsim) {
  
  print(paste("simulation",i, sep = " "))
  
  set.seed(i^3+5444)
  Z = matrix(runif(N * p_Clin, 0, 1), nrow = N, ncol = p_Clin)
  Z[,1] <- sample(0:3, N, replace = TRUE)
  Z[,2] <- sample(0:1, N, replace = TRUE)
  colnames(Z) = paste0("z", seq(1, p_Clin))
  Z <- data.frame(Z)
  
  X = mvtnorm::rmvnorm(N, sigma = CorrX)
  colnames(X) = paste0("x", seq(1, p))
  
  prob <- logit_prob(z = Z, x = X, betas = betas)
  Y <- rbinom(N, size = 1, prob = prob)
  
  #specify folds for hyperparameter tuning
  folds <- CVfolds(Y=Y, model = "logistic", balance = F, kfold = 5, nrepeat = 1)
  
  ### 1. Random Forest (classification) ###
  DFtrain <- data.frame(Ydf = factor(Y), Xdf = cbind(X, Z))
  DFtest  <- data.frame(Ydf = factor(Ytest), Xdf = cbind(Xtest, Ztest))
  
  RF <- rfsrc(Ydf ~ ., data = DFtrain, ntree = 500, mtry = floor((p + p_Clin) / 3),
              var.used = "all.trees", importance = "none", family = "class")
  
  preds_RF_prob <- predict.rfsrc(RF, newdata = DFtest)$predicted[, 2]  # probability for class 1
  AUC_RF[i] <- auc(Ytest, preds_RF_prob)
  Brier_RF[i] <- mean((preds_RF_prob - Ytest)^2)
  remove(RF, preds_RF_prob, DFtrain, DFtest)
  
  ### 2. Gradient Boosting (XGBoost) ###
  
  tune_grid <- expand.grid(
    nrounds = c(100, 200),
    eta = c(0.01, 0.1),
    max_depth = c(2,4),
    gamma = 0,
    colsample_bytree = 0.75,
    min_child_weight = N/15,
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
    x = cbind(X, Z),
    y = Y_class,
    trControl = tune_control,
    tuneGrid = tune_grid,
    method = "xgbTree",
    verbosity = 0,
    metric = "ROC"
  )
  xgb_tune$bestTune
  
  xgb_pred_prob <- predict(xgb_tune, cbind(Xtest, Ztest), type = "prob")[, "yes"]
  AUC_GB[i] <- auc(Ytest, xgb_pred_prob)
  Brier_GB[i] <- mean((xgb_pred_prob - Ytest)^2)
  remove(xgb_tune, xgb_pred_prob)
  
  ### 3. Ridge Regression (binomial, clinical unpenalized) ###
  Lam1 <- porridge::optPenaltyGLM.kCVauto(
    Y = Y, X = X, U = as.matrix(Z),
    lambdaInit = 10, lambdaGinit = 0,
    Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
    model = "logistic", folds = folds)
  
  if (Lam1 == Inf) Lam1 <- 1e8
  
  RidgeFit <- ridgeGLM(
    Y = Y, U = as.matrix(Z), X = X,
    lambda = Lam1, lambdaG = 0,
    Dg = matrix(0, ncol = ncol(X), nrow = ncol(X)),
    model = "logistic"
  )
  
  ridge_prob <- (cbind(as.matrix(Ztest), Xtest) %*% RidgeFit)[, 1]
  ridge_prob <- 1 / (1 + exp(-ridge_prob))  # convert logit to probability
  AUC_ridge[i] <- auc(Ytest, ridge_prob)
  Brier_ridge[i] <- mean((ridge_prob - Ytest)^2)
  remove(RidgeFit, ridge_prob)
  
  ### 4. Lasso (logistic) ###
  Las <- cv.glmnet(cbind(as.matrix(Z), X), Y, alpha = 1, nfolds = 10, family = "binomial",
                   penalty.factor = c(rep(0, p_Clin), rep(1, p)))$lambda.min
  
  LassoFit <- glmnet(cbind(as.matrix(Z), X), Y, alpha = 1, lambda = Las, family = "binomial",
                     penalty.factor = c(rep(0, p_Clin), rep(1, p)))
  
  lasso_prob <- predict(LassoFit, newx = cbind(as.matrix(Ztest), Xtest),
                        s = Las, type = "response")[, 1]
  AUC_Lasso[i] <- auc(Ytest, lasso_prob)
  Brier_Lasso[i] <- mean((lasso_prob - Ytest)^2)
  remove(LassoFit, lasso_prob)
  
  
  ####### 6. Fit fusedTree ######
  dat=cbind.data.frame(Y,Z)
  # fitting tree
  rp <- rpart(Y~.,data = dat, control = rpart.control(xval = 5,  maxdepth = 2),
              model = T)
  cp = rp$cptable[,1][which.min(rp$cptable[,4])]
  Treefit <- rp
  remove(dat,rp)
  n_leaves_vec[i] <- sum(Treefit$frame$var == "<leaf>")
  
  ####### 4. Fit fully FusedReg Tree ######
  
  ## FusedTree ##
  tryCatch({
    
    foldsTree <- CVfoldsTree(Y = Y, Tree = Treefit, Z = data.frame(Z), model = "logistic", nrepeat = 1, kfold = 5)
    
    # 6a. Regular fusedTree
    optPenalties <- PenOpt(Treefit, X, Y, data.frame(Z), "logistic", lambdaInit = 10, alphaInit = 10,
                           loss = "loglik", LinVars = FALSE, folds = foldsTree, multistart = FALSE)
    
    fit <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "logistic",
                     lambda = optPenalties[1], alpha = optPenalties[2], verbose = FALSE)
    
    Ypred_prob <- predict(fit, newX = Xtest, newY = Ytest, newZ = Ztest)$Ypred
    AUC_fusedTree[i] <- as.numeric(pROC::auc(Ytest, Ypred_prob))
    Brier_fusedTree[i] <- mean((Ytest - Ypred_prob)^2)
    alphas[i] <- optPenalties[2]
    remove(fit, Ypred_prob, optPenalties)
    
    # 6b. Zero fused
    optPenalties <- PenOpt(Treefit, X, Y, data.frame(Z), "logistic", lambdaInit = 10, alphaInit = 0,
                           loss = "loglik", LinVars = FALSE, folds = foldsTree, multistart = FALSE)
    
    fit1 <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "logistic",
                      lambda = optPenalties[1], alpha = 0, verbose = F)
    Ypred_prob <- predict(fit1, newX = Xtest, newY = Ytest, newZ = Ztest)$Ypred
    AUC_zeroFus[i] <- as.numeric(pROC::auc(Ytest, Ypred_prob))
    Brier_zeroFus[i] <- mean((Ytest - Ypred_prob)^2)
    remove(fit1, optPenalties, Ypred_prob)
    
    # 6c. Full fused (large alpha)
    optPenalties <- PenOpt(Treefit, X, Y, data.frame(Z), "logistic", lambdaInit = 10, alphaInit = 10e10,
                           loss = "loglik", LinVars = FALSE, folds = foldsTree, multistart = FALSE)
    
    fit2 <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "logistic",
                      lambda = optPenalties[1], alpha = 10e10, verbose = F)
    Ypred_prob <- predict(fit2, newX = Xtest, newY = Ytest, newZ = Ztest)$Ypred
    AUC_fulFus[i] <- as.numeric(pROC::auc(Ytest, Ypred_prob))
    Brier_fulFus[i] <- mean((Ytest - Ypred_prob)^2)
    remove(fit2, optPenalties, Ypred_prob)
    
  }, error = function(e) {
    warning(paste("FusedTree error in simulation", i, ":", e$message))
    
    # fallback: use ridge probabilities already computed above
    AUC_fusedTree[i] <- NA
    Brier_fusedTree[i] <- NA
    
    AUC_zeroFus[i] <- NA
    Brier_zeroFus[i] <- NA
    
    AUC_fulFus[i] <- NA
    Brier_fulFus[i] <- NA
    
    alphas[i] <- NA
  })
  
  ####### 5. Fit glinternet ######
  numLevs <- c(4, 2, 1, 1, 1, rep(1, ncol(X)))  # set factor levels for categorical vars
  
  result <- tryCatch({
    Lam <- glinternet.cv(
      X = cbind(as.matrix(Z), X), 
      Y = Y, 
      numLevels = numLevs, 
      nFolds = 5,
      interactionCandidates = 1:5,
      family = "binomial"
    )
    
    fit <- glinternet(
        X = cbind(as.matrix(Z), X), 
        Y = Y, 
        numLevels = numLevs, 
        lambda = Lam$lambdaHat,
        interactionCandidates = 1:5,
        family = "binomial"
    )
    
    Ypred_prob <- predict(fit, X = cbind(as.matrix(Ztest), Xtest), 
                          type = "response", lambda = Lam$lambdaHat)[,1]
    
    AUC_glinternet[i] <- as.numeric(pROC::auc(Ytest, Ypred_prob))
    Brier_glinternet[i] <- mean((Ytest - Ypred_prob)^2)
    
  }, error = function(e) {
    message(sprintf("Simulation %d: glinternet failed or timed out", i))
    AUC_glinternet[i] <- NA
    Brier_glinternet[i] <- NA
  })
  
  # Cleanup
  remove(Lam, fit, Ypred_prob)
  
  remove(X,Y,Z,Treefit,folds, foldsTree)
  gc()
}

Brier_fulFus[Brier_fulFus == 0] <- NA
Brier_fusedTree[Brier_fusedTree == 0] <- NA
Brier_zeroFus[Brier_zeroFus == 0] <- NA
AUC_fulFus[AUC_fulFus == 0] <- NA
AUC_fusedTree[AUC_fusedTree == 0] <- NA
AUC_zeroFus[AUC_zeroFus == 0] <- NA

AUC_df <- data.frame(
  fusedTree  = AUC_fusedTree,
  fullFus    = AUC_fulFus,
  ZeroFus    = AUC_zeroFus,
  GB         = AUC_GB,
  RF         = AUC_RF,
  Ridge      = AUC_ridge,
  Lasso      = AUC_Lasso,
  Glinternet = AUC_glinternet
)

Brier_df <- data.frame(
  fusedTree  = Brier_fusedTree,
  fullFus    = Brier_fulFus,
  ZeroFus    = Brier_zeroFus,
  GB         = Brier_GB,
  RF         = Brier_RF,
  Ridge      = Brier_ridge,
  Lasso      = Brier_Lasso,
  Glinternet = Brier_glinternet
)



colMeans(AUC_df, na.rm = T)
colMeans(Brier_df,na.rm = T)
results <- list(AUC = AUC_df, Brier = Brier_df)
nm = paste(N,Nsim,"Binary_Interaction.Rdata",sep = "_")

save(results,file = nm)



getwd()

## do glinternet later as something is wrong
N <- 300 # N = 100 or N = 200 do not work unfortunately
AUC_glinternet <- numeric(Nsim)
Brier_glinternet <- numeric(Nsim)
library(R.utils)


for (i in 1:Nsim) {
  
  print(paste("simulation", i, sep = " "))
  
  set.seed(i^3+5444)
  Z = matrix(runif(N * p_Clin, 0, 1), nrow = N, ncol = p_Clin)
  Z[,1] <- sample(0:3, N, replace = TRUE)
  Z[,2] <- sample(0:1, N, replace = TRUE)
  colnames(Z) = paste0("z", seq(1, p_Clin))
  Z <- data.frame(Z)
  
  X = mvtnorm::rmvnorm(N, sigma = CorrX)
  colnames(X) = paste0("x", seq(1, p))
  
  prob <- logit_prob(z = Z, x = X, betas = betas)
  Y <- rbinom(N, size = 1, prob = prob)
  
  ####### 5. Fit glinternet ######
  numLevs <- c(4, 2, 1, 1, 1, rep(1, ncol(X)))  # set factor levels for categorical vars
  
  result <- tryCatch({
    Lam <- glinternet.cv(
      X = cbind(as.matrix(Z), X), 
      Y = Y, 
      numLevels = numLevs, 
      nFolds = 5,
      interactionCandidates = 1:5,
      family = "binomial"
    )
    
    fit <- withTimeout({
      glinternet(
        X = cbind(as.matrix(Z), X), 
        Y = Y, 
        numLevels = numLevs, 
        lambda = Lam$lambdaHat,
        interactionCandidates = 1:5,
        family = "binomial"
      )
    }, timeout = 60, onTimeout = "error")
    
    Ypred_prob <- predict(fit, X = cbind(as.matrix(Ztest), Xtest), 
                          type = "response", lambda = Lam$lambdaHat)[,1]
    
    AUC_glinternet[i] <- as.numeric(pROC::auc(Ytest, Ypred_prob))
    Brier_glinternet[i] <- mean((Ytest - Ypred_prob)^2)
    
  }, error = function(e) {
    message(sprintf("Simulation %d: glinternet failed or timed out", i))
    AUC_glinternet[i] <- NA
    Brier_glinternet[i] <- NA
  })
  
  # Cleanup
  remove(Lam, fit, Ypred_prob)
  remove(X, Y, Z)
  gc()
}
load("300_500_Binary_Interaction.Rdata")
AUC_df <- results[[1]]
Brier_df <- results[[2]]

AUC_df <- cbind.data.frame(AUC_df, Glinternet = AUC_glinternet)
Brier_df <- cbind.data.frame(Brier_df, Glinternet = Brier_glinternet)

colMeans(AUC_df, na.rm = T)
colMeans(Brier_df,na.rm = T)
results <- list(AUC = AUC_df, Brier = Brier_df)
nm = paste(N,Nsim,"Binary_Interaction.Rdata",sep = "_")
save(results,file = nm)

library(tidyverse)
load("100_500_Binary_Interaction.Rdata")
# Add identifier columns
AUC_df_100 <- results$AUC
Brier_df_100 <- results$Brier
colMeans(AUC_df_100, na.rm = T)
colMeans(Brier_df_100, na.rm = T)
AUC_100_long <- pivot_longer(AUC_df_100, everything(), names_to = "Model", values_to = "Score") %>%
  mutate(SampleSize = "N = 100", Metric = "AUC")

Brier_100_long <- pivot_longer(Brier_df_100, everything(), names_to = "Model", values_to = "Score") %>%
  mutate(SampleSize = "N = 100", Metric = "Brier")

load("300_500_Binary_Interaction.Rdata")
AUC_df_300 <- results$AUC
Brier_df_300 <- results$Brier
colMeans(AUC_df_300, na.rm = T)
colMeans(Brier_df_300, na.rm = T)

AUC_300_long <- pivot_longer(AUC_df_300, everything(), names_to = "Model", values_to = "Score") %>%
  mutate(SampleSize = "N = 300", Metric = "AUC")

Brier_300_long <- pivot_longer(Brier_df_300, everything(), names_to = "Model", values_to = "Score") %>%
  mutate(SampleSize = "N = 300", Metric = "Brier")

# Combine all
perf_df <- bind_rows(AUC_100_long, Brier_100_long, AUC_300_long, Brier_300_long)

perf_df <- perf_df %>%
  mutate(
    Model = recode(Model,
                   fusedTree = "FusedTree",
                   fullFus = "FullFus"),  # Rename here
    Model = factor(Model, levels = c("FusedTree", "FullFus", "ZeroFus", 
                                     "GB", "RF", "Ridge", "Lasso", "Glinternet"))
  )
# omit Glinternet results for N=100 because these are zero
perf_df <- perf_df %>% filter(!(Model == "Glinternet" & SampleSize == "N = 100"))
## plot

library(viridis)
name <- paste("Figure_Simulation_Binary.pdf")
pdf(name, width=8,height=5)
ggplot(perf_df, aes(x = Model, y = Score, fill = Model)) +
  geom_boxplot(outlier.shape = NA,fatten=2) +
  facet_grid(Metric ~ SampleSize, scales = "free_y") +
  theme_light(base_size = 14) + scale_fill_viridis(discrete = T, option = "D") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       y = "Score",
       x = "Model")
dev.off()


