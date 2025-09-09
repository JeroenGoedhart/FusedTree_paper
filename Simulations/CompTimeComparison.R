setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
gc()
library(VGAM)
library(glmnet)
library(randomForestSRC)
library(penalized)
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
  1*(z[,1] < 1.5)*(z[,2] < 0.5)*(-10 + x[,1:(p/4)] %*% betas[1:(p/4)]*8) + #node1
    1*(z[,1] < 1.5)*(z[,2] >= 0.5)*(-2 + x[,1:(p/4)]%*% betas[1:(p/4)]*2) + #node2
    1*(z[,1] >= 1.5)*(z[,3] < 0.5)*(5 + x[,1:(p/4)]%*% betas[1:(p/4)]*0.5) +
    1*(z[,1] >= 1.5)*(z[,3] >= 0.5)*(10 + x[,1:(p/4)]%*% betas[1:(p/4)]*0.125) +
    x[,(p/4+1):p]%*% betas[(p/4+1):p] + 5*z[,4]
}

N_values <- c(100, 200, 500)
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
    
    colnames(time_mat) <- c("Ridge", "FusedTree", "RF", "GB", "Lasso", "Glinternet")
    set.seed(p)
    ids <- sample(1:ncol(CorrX), p, replace = FALSE)
    CorrX1 <- CorrX[ids, ids]
    
    for (i in 1:Nsim) {
      
      print(paste("Running repeat: ", i))
      
      set.seed(i ^ 3 + 3)
      Z <- matrix(runif(N * 5, 0, 1), nrow = N)
      Z <- data.frame(Z)
      colnames(Z) <- paste0("z", seq_len(ncol(Z)))
      Z[,1] <- as.integer(sample(0:3, N, replace = TRUE))
      Z[,2] <- as.integer(sample(0:1, N, replace = TRUE))
      
      X <- mvtnorm::rmvnorm(N, sigma = CorrX1)
      colnames(X) <- paste0("x", seq_len(ncol(X)))
      
      betas <- matrix(rlaplace(p, 0, 20 / p), ncol = 1)
      Y <- g(z = Z, x = X, betas = betas) + rnorm(N, 0, 1)
      Y <- Y[,1]
      Z$z1 <- as.factor(Z$z1)
      
      folds <- multiridge::CVfolds(Y = Y, model = "linear", 
                                   balance = FALSE, kfold = 5, nrepeat = 1)
      
      ### Ridge ###
      tryCatch({
        time_mat[i, "Ridge"] <- suppressMessages(suppressWarnings(system.time({
          
          lambda_opt <- optL2(response = Y, penalized = X, 
                              unpenalized = model.matrix(~., Z)[,-1],
                              trace = FALSE, fold = 5)
          fit <- penalized(response = Y, penalized = X, 
                           unpenalized = model.matrix(~., Z)[,-1],
                           lambda1 = 0, lambda2 = 10$lambda,
                           trace = FALSE)
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("Ridge error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "Ridge"] <- NA
      })
      
      ### FusedTree ###
      tryCatch({
        time_mat[i, "FusedTree"] <- suppressMessages(suppressWarnings(system.time({
          
          dat <- cbind.data.frame(Y, data.frame(Z))
          Treefit <- rpart(Y ~ ., data = dat, 
                           control = rpart.control(xval = 5, maxdepth = 3), 
                           model = TRUE)
          
          n_leaves_vec[i] <- sum(Treefit$frame$var == "<leaf>")
          
          foldsTree <- CVfoldsTree(Y = Y, Tree = Treefit, Z = data.frame(Z), 
                                   model = "linear", nrepeat = 1, kfold = 5)
          
          optPenalties <- PenOpt(Treefit, X, Y, data.frame(Z), "linear", 
                                 lambdaInit = 10, alphaInit = 10,
                                 loss = "loglik", LinVars = FALSE, 
                                 folds = foldsTree, multistart = FALSE)
          
          fit <- fusedTree(Treefit, X, Y, data.frame(Z), LinVars = FALSE, 
                           model = "linear",
                           lambda = optPenalties[1], alpha = optPenalties[2])
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("fusedTree error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "FusedTree"] <- NA
      })
      
      ### Random Forest ###
      tryCatch({
        time_mat[i, "RF"] <- suppressMessages(suppressWarnings(system.time({
          
          DFtrain <- data.frame(Ydf = Y, Xdf = cbind(X, Z))
          
          sink(tempfile())  # redirect printed output
          invisible(capture.output({
            tune_result <- tune.nodesize(Ydf ~ ., data = DFtrain,
                                         nodesizeTry = c(3, 5, 10, 50),
                                         ntreeTry = 500,
                                         doBest = FALSE)
          }))
          sink()
          
          RF <- rfsrc(Ydf ~ ., data = DFtrain,ntree = 500, mtry = p / 3,
                      nodesize = tune_result$nsize.opt)
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("RF error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "RF"] <- NA
      })
      
      ### GB ###
      tryCatch({
        time_mat[i, "GB"] <- suppressMessages(suppressWarnings(system.time({
          
          tune_grid <- expand.grid(nrounds = c(200, 500), eta = c(0.01 ,0.1),
                                   max_depth = 3, gamma = 0, colsample_bytree = .75,
                                   min_child_weight = 30, subsample = 1)
          
          tune_control <- caret::trainControl(method = "cv", number = 5, verboseIter = FALSE)
          
          xgb_tune <- caret::train(x = cbind(model.matrix( ~ .,Z), X), y = Y,
                                   trControl = tune_control,
                                   tuneGrid = tune_grid, method = "xgbTree",
                                   verbose = FALSE, verbosity = 0)
        })["elapsed"]))
      }, error = function(e) {
        warning(sprintf("GB error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "GB"] <- NA
      })
      
      ### Lasso ###
      tryCatch({
        time_mat[i, "Lasso"] <- system.time({
          
          penfact <- c(rep(0, ncol(model.matrix(~., Z))), rep(1, ncol(X)))
          
          lam <- cv.glmnet(cbind(model.matrix( ~ ., Z), X), Y, alpha = 1, nfolds = 10,
                           penalty.factor = penfact)$lambda.min
          
          LassoFit <- glmnet(cbind(model.matrix( ~ ., Z),X), Y, alpha = 1, 
                             lambda = lam, penalty.factor = penfact)
          
        })["elapsed"]
      }, error = function(e) {
        warning(sprintf("Lasso error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "Lasso"] <- NA
      })
      
      ### glinternet ###
      tryCatch({
        time_mat[i, "Glinternet"] <- system.time({
          
          numLevs <- c(4, 2, 1, 1, 1, rep(1, ncol(X)))
          Z$z1 <- as.numeric(Z$z1) - 1
          
          Lam <- glinternet.cv(X = cbind(as.matrix(Z),X), Y = Y, numLevels = numLevs,
                               nFolds = 5, interactionCandidates = c(1:5),
                               family = "gaussian")
          
          fit <- glinternet(X = cbind(as.matrix(Z),X), Y = Y, numLevels = numLevs,
                            lambda = Lam$lambdaHat,
                            interactionCandidates = c(1:5),
                            family = "gaussian")
        })["elapsed"]
      }, error = function(e) {
        warning(sprintf("glinternet error at sim %d (N=%d, p=%d): %s", i, N, p, conditionMessage(e)))
        time_mat[i, "Glinternet"] <- NA
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
      summarise(median_time = median(time, na.rm = TRUE), .groups = "drop") %>%
      mutate(N = N, p = p)
    
    leaf_df <- data.frame(sim = 1:Nsim, leaves = n_leaves_vec, N = N, p = p)
    
    results <- rbind(results, time_df)
    
    results_median <- rbind(results_median, time_df1)
    
    res_leafs <- rbind(res_leafs, leaf_df)
  }
}
warnings()
allres1 <- list(time_average = results,time_median = results_median, leaves = res_leafs)
save(allres1, file = "CompTime_Results_New.Rdata")

##### Plotting #####

load("CompTime_Results.Rdata")
results <- allres$time
library(ggplot2); library(ggsci)

name <- "CompTimeComparison.pdf"
pdf(name, width=7,height=4)
ggplot(results, aes(x = p, y = mean_time, color = model, fill = model, shape = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2, stroke = 0.8) +
  facet_wrap(
    ~ N,
    scales = "free_y",
    labeller = as_labeller(function(x) paste0("N = ", x)),
    nrow = 1,
    ncol = 3
  ) +
  theme_light()  +
  theme(legend.position = "bottom") + 
  scale_color_jama() +
  scale_fill_jama() +
  scale_shape_manual(values = c(21,4, 22:25)) +
  scale_x_continuous(trans = 'log10') +
  labs(
    x = "Number of omics variables (p)",
    y = "Average computation time (s)",
    title = "",
    color = "Model",
    fill = "Model",
    shape = "Model"
  )
dev.off()
