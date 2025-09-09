setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
library(mvtnorm)
library(Matrix)
library(VGAM)
library(rpart); library(rpart.plot)


N <- 300; p <- 500; p_Clin <- 5
Ntest <- 5000
set.seed(4)

g <- function (z,x,betas){
  1*(z[,1] < 1.5)*(z[,2] < 0.5)*(-10 + x[,1:(p/4)] %*% betas[1:(p/4)]*0) + #node1
    1*(z[,1] < 1.5)*(z[,2] >= 0.5)*(-5 + x[,1:(p/4)]%*% betas[1:(p/4)]*2.5) + #node2
    1*(z[,1] >= 1.5)*(z[,3] < 0.5)*(5 + x[,1:(p/4)]%*% betas[1:(p/4)]*3) +
    1*(z[,1] >= 1.5)*(z[,3] >= 0.5)*(10 + x[,1:(p/4)]%*% betas[1:(p/4)]*4) +
    x[,(p/4+1):p] %*% betas[(p/4+1):p]
}

set.seed(344)
betas <- matrix(rlaplace(p, 0, 30/p), ncol = 1)
Xtest = mvtnorm::rmvnorm(Ntest, sigma = diag(1,p))
colnames(Xtest) = paste0("x",seq(1,p))
Ztest = matrix(runif(Ntest*p_Clin,0,1),nrow = Ntest, ncol = p_Clin)
Ztest[,1] <- sample(0:3,Ntest,replace = T)
Ztest[,2] <- sample(0:1,Ntest,replace = T)
Ztest = data.frame(Ztest)
colnames(Ztest) = paste0("z",seq(1,p_Clin))
set.seed(4)
Ytest = g(z = Ztest, x = Xtest, betas=betas)[,1] + rnorm(Ntest,0,1)
var(Ytest)
Ztest$z1 <- as.factor(Ztest$z1)
Ztest$z2 <- as.factor(Ztest$z2)

MSE_symFus <- c()
MSE_asymFus <- c()
Nsim <- 500
for (i in 1:Nsim) {
  
  print(paste("simulation",i, sep = " "))
  ## simulating clinical covariates
  set.seed(i^3+1)
  Z <- matrix(runif(N*5, 0, 1), nrow = N)
  Z[,1] <- sample(0:3, N, replace = TRUE)
  Z[,2] <- sample(0:1, N, replace = TRUE)
  colnames(Z) <- paste0("z", seq_len(ncol(Z)))
  Z <- data.frame(Z)
  X <- mvtnorm::rmvnorm(N, sigma = diag(1,p))
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  Y <- g(z = Z, x = X, betas = betas) + rnorm(N, 0, 1)
  Y <- Y[,1]
  Z$z1 <- as.factor(Z$z1)
  Z$z2 <- as.factor(Z$z2)
  dat <- cbind.data.frame(Y, data.frame(Z))
  rp <- rpart(Y ~ ., data = dat, control = rpart.control(xval = 5, minbucket = 20), model = TRUE)
  cp <- rp$cptable[,1][which.min(rp$cptable[,4])]
  Treefit <- prune(rp, cp = cp)
  remove(dat,rp,cp)
  foldsTree <- CVfoldsTree(Y = Y, Tree = Treefit, Z = data.frame(Z), model = "linear", nrepeat = 1, kfold = 5)
  
  ## asymmetric fusion ##
  tryCatch({
    optPenalties <- PenOpt(Tree = Treefit, X = X, Y = Y, Z = data.frame(Z), "linear", 
                           lambdaInit = 100, alphaInit = 1000, multistart = FALSE,
                           loss = "loglik", LinVars = FALSE, folds = foldsTree,
                           symFusion = FALSE)
    
    fit <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "linear",
                     lambda = optPenalties[1], alpha = optPenalties[2], symFusion = FALSE)
    
    Preds <- predict.fusedTree(fit, newX = Xtest, newZ = Ztest, newY = Ytest)
    MSE_asymFus[i] <- mean((Preds$Ypred - Ytest)^2)
    
    remove(optPenalties, fit, Preds)
  }, error = function(e) {
    warning(paste("Asymmetric fusion failed at iteration", i, ":", e$message))
    MSE_asymFus[i] <- NA
  })
  
  ## symmetric fusion ##
  tryCatch({
    optPenalties1 <- PenOpt(Tree = Treefit, X = X, Y = Y, Z = data.frame(Z), "linear", 
                            lambdaInit = 100, alphaInit = 1000, multistart = TRUE,
                            loss = "loglik", LinVars = FALSE, folds = foldsTree,
                            symFusion = TRUE)
    
    fit1 <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "linear",
                      lambda = optPenalties1[1], alpha = optPenalties1[2], symFusion = TRUE)
    
    Preds1 <- predict.fusedTree(fit1, newX = Xtest, newZ = Ztest, newY = Ytest)
    MSE_symFus[i] <- mean((Preds1$Ypred - Ytest)^2)
    
    remove(optPenalties1, fit1, Preds1)
  }, error = function(e) {
    warning(paste("Symmetric fusion failed at iteration", i, ":", e$message))
    MSE_symFus[i] <- NA
  })
  
  remove(optPenalties1,fit1,Preds1,cp,Treefit,rp)
  remove(X,Z,Y,foldsTree)
  gc()
}

mean(MSE_asymFus, na.rm = T)
mean(MSE_symFus, na.rm = T)

mean(MSE_asymFus)/mean(MSE_symFus)


fit$Effects[1:9]


betas[1]*c(8,2,0.5,0.125)
View(cbind(fit$Effects[5:8],fit1$Effects[5:8]))
   
  