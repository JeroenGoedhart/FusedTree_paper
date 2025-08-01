setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
gc()
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
library(dplyr)
library(tidyr)
library(microbenchmark)  # optional, but more precise timing


load("correlationmatrix.Rdata")

g <- function (z,x,betas){
  1*(z[,1] < 1.5)*(z[,2] < 0.5)*(-10 + x[,1:p/4] %*% betas[1:p/4]*8) + #node1
    1*(z[,1] < 1.5)*(z[,2] >= 0.5)*(-2 + x[,1:p/4]%*% betas[1:p/4]*2) + #node2
    1*(z[,1] >= 1.5)*(z[,3] < 0.5)*(5 + x[,1:p/4]%*% betas[1:p/4]*0.5) +
    1*(z[,1] >= 1.5)*(z[,3] >= 0.5)*(10 + x[,1:p/4]%*% betas[1:p/4]*0.125) +
    x[,(p/4+1):p]%*% betas[(p/4+1):p] + 5*z[,4]
}
N_values <- c(100,200)
p_values <- c(10,20)
Nsim <- 10

N_values <- c(100, 200)
p_values <- c(50, 100, 200, 500, 1000)
Nsim <- 50

results <- data.frame()
results_median <- data.frame()
res_leafs <- data.frame()
for (N in N_values) {
  for (p in p_values) {
    message(sprintf("Running N = %d, p = %d", N, p))
    
    time_mat <- matrix(0, nrow = Nsim, ncol = 6)
    n_leaves_vec <- c()  # to store number of leaves per simulation
    
    colnames(time_mat) <- c("ridge", "fusedTree", "rf", "gb", "lasso", "glinternet")
    set.seed(p)
    ids <- sample(1:ncol(CorrX), p, replace = FALSE)
    CorrX1 <- CorrX[ids, ids]
    
    for (i in 1:Nsim) {
      set.seed(i^2 + 3)
      Z <- matrix(runif(N*5, 0, 1), nrow = N)
      Z[,1] <- sample(0:3, N, replace = TRUE)
      Z[,2] <- sample(0:1, N, replace = TRUE)
      
      colnames(Z) <- paste0("z", seq_len(ncol(Z)))
      Z <- data.frame(Z)
      
      X <- mvtnorm::rmvnorm(N, sigma = CorrX1)
      colnames(X) <- paste0("x", seq_len(ncol(X)))
      
      betas <- matrix(rlaplace(p, 0, 10/p), ncol = 1)
      Y <- g(z = Z, x = X, betas = betas) + rnorm(N, 0, 1)
      Y <- Y[,1]
      Z$z1 <- as.factor(Z$z1)
      folds <- multiridge::CVfolds(Y = Y, model = "linear", balance = FALSE, kfold = 5, nrepeat = 1)
      
      ### Ridge ###
      tryCatch({
        time_mat[i, "ridge"] <- suppressMessages(suppressWarnings(system.time({
          Lam1 <- porridge::optPenaltyGLM.kCVauto(Y = Y, X = X, U = model.matrix(~.,Z),
                                                  lambdaInit = 10, lambdaGinit = 0,
                                                  Dg = matrix(0, ncol=ncol(X), nrow=ncol(X)),
                                                  model = "linear", folds = folds, loss = "sos")
          if (Lam1 == Inf) Lam1 <- 10e8
          fit <- ridgeGLM(Y = Y, U = model.matrix(~.,Z), X = X, lambda = Lam1, lambdaG = 0,
                          Dg = matrix(0, ncol=ncol(X), nrow=ncol(X)), model = "linear")
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("Ridge error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "ridge"] <- NA
      })
      
      ### FusedTree ###
      tryCatch({
        time_mat[i, "fusedTree"] <- suppressMessages(suppressWarnings(system.time({
          dat <- cbind.data.frame(Y, data.frame(Z))
          Treefit <- rpart(Y ~ ., data = dat, control = rpart.control(xval = 5, maxdepth = 3), model = TRUE)
          n_leaves_vec[i] <- sum(Treefit$frame$var == "<leaf>")
          foldsTree <- CVfoldsTree(Y = Y, Tree = Treefit, Z = data.frame(Z), model = "linear", nrepeat = 1, kfold = 5)
          optPenalties <- PenOpt(Treefit, X, Y, data.frame(Z), "linear", lambdaInit = 10, alphaInit = 10,
                                 loss = "loglik", LinVars = FALSE, folds = foldsTree, multistart = FALSE)
          fit <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, model = "linear",
                           lambda = optPenalties[1], alpha = optPenalties[2])
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("fusedTree error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "fusedTree"] <- NA
      })
      
      ### Random Forest ###
      tryCatch({
        time_mat[i, "rf"] <- suppressMessages(suppressWarnings(system.time({
          DFtrain <- data.frame(Ydf = Y, Xdf = cbind(X, Z))
          sink(tempfile())  # redirect printed output
          invisible(capture.output({
            tune_result <- tune.nodesize(Ydf ~ ., data = DFtrain,
                                         nodesizeTry = c(3, 5, 10, 50),
                                         ntreeTry = 500,
                                         doBest = FALSE)
          }))
          sink()  # 
          
          RF <- rfsrc(Ydf ~ .,data=DFtrain,ntree = 500, mtry = p/4,
                      nodesize = tune_result$nsize.opt)
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("RF error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "rf"] <- NA
      })
      
      ### GB ###
      tryCatch({
        time_mat[i, "gb"] <- suppressMessages(suppressWarnings(system.time({
          tune_grid <- expand.grid(nrounds = c(200, 500), eta = c(0.01,0.1), max_depth = 3,
                                   gamma = 0, colsample_bytree = 0.75,
                                   min_child_weight = 30, subsample = 1)
          tune_control <- caret::trainControl(method = "cv", number = 5, verboseIter = FALSE)
          xgb_tune <- caret::train(x = cbind(model.matrix(~.,Z),X), y = Y,
                                   trControl = tune_control,
                                   tuneGrid = tune_grid, method = "xgbTree",
                                   verbose = FALSE, verbosity = 0)
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("GB error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "gb"] <- NA
      })
      
      ### Lasso ###
      tryCatch({
        time_mat[i, "lasso"] <- system.time({
          
          lam <- cv.glmnet(cbind(model.matrix(~., Z),X), Y, alpha = 1, nfolds = 10,
                           penalty.factor = c(rep(0, ncol(model.matrix(~., Z))), rep(1, ncol(X))))$lambda.min
          LassoFit <- glmnet(cbind(model.matrix(~., Z),X), Y, alpha = 1, lambda = lam,
                             penalty.factor = c(rep(0, ncol(model.matrix(~., Z))), rep(1, ncol(X))))
          LassoFit$beta
        })["elapsed"]
      }, error = function(e) {
        warning(sprintf("Lasso error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "lasso"] <- NA
      })
      
      ### glinternet ###
      tryCatch({
        time_mat[i, "glinternet"] <- system.time({
          numLevs <- c(4, 2, 1, 1, 1, rep(1, ncol(X)))
          Z$z1 <- as.numeric(Z$z1) - 1
          Lam <- glinternet.cv(X = cbind(as.matrix(Z),X), Y = Y, numLevels = numLevs,
                               nFolds = 5, interactionCandidates = 1:5,
                               family = "gaussian")
          fit <- glinternet(X = cbind(as.matrix(Z),X), Y = Y, numLevels = numLevs,
                            lambda = Lam$lambdaHat,
                            interactionCandidates = 1:5,
                            family = "gaussian")
        })["elapsed"]
      }, error = function(e) {
        warning(sprintf("glinternet error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "glinternet"] <- NA
      })
    }
    time_df <- as.data.frame(time_mat) %>%
      mutate(sim = 1:Nsim) %>%
      pivot_longer(cols = -sim, names_to = "model", values_to = "time") %>%
      group_by(model) %>%
      summarise(mean_time = mean(time, na.rm = TRUE), .groups = "drop") %>%
      mutate(N = N, p = p)
    
    time_df1 <- as.data.frame(time_mat) %>%
      mutate(sim = 1:Nsim) %>%
      pivot_longer(cols = -sim, names_to = "model", values_to = "time") %>%
      group_by(model) %>%
      summarise(median_time = mean(time, na.rm = TRUE), .groups = "drop") %>%
      mutate(N = N, p = p)
    
    leaf_df <- data.frame(sim = 1:Nsim, leaves = n_leaves_vec, N = N, p = p)
    
    results <- rbind(results, time_df)
    results_median <- rbind(results_median, time_df1)
    res_leafs <- rbind(res_leafs, leaf_df)
  }
}

warnings()
allres <- list(time = results, leaves = res_leafs)
save(allres, file = "comptimeres.Rdata")


library(ggplot2); library(ggsci)

name <- "CompTimeComparison.pdf"
pdf(name, width=7,height=3)
ggplot(results, aes(x = p, y = mean_time, color = model)) +
  geom_line(linewidth= 1.2) +
  geom_point() +
  facet_wrap(~N, scales = "free_y", labeller = as_labeller(function(x) paste0("N = ", x))) +
  theme_light() + scale_color_jama() + 
  scale_x_continuous(trans = 'log10') +
  labs(x = "Number of omics variables (p)",
       y = "Average computation time (s)",
       title = "",
       color = "Model")
dev.off()


subset_results <- subset(results, model %in% c("lasso", "ridge"))

# Plot mean_time vs p, colored by model, faceted by N
ggplot(subset_results, aes(x = p, y = mean_time, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_point() +
  facet_wrap(~ N, scales = "free") +
  labs(title = "Computation Time for Lasso and Ridge",
       x = "Number of predictors (p)",
       y = "Mean Time (s)",
       color = "Model") +
  theme_minimal()
