setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
gc()
library(partykit)
library(mvtnorm)
library(VGAM)  # for rlaplace
library(MASS)  # for predict.lm fallback in mob
library(dplyr)
library(fusedTree)
library(rpart)
library(rpart.plot)



# Define custom data generation function
g <- function(z, x, betas) {
  p <- ncol(x)
  1 * (z[,1] < 1.5) * (z[,2] < 0.5) * (-10 + x[,1:(p/2)] %*% betas[1:(p/2)] * 8) +
    1 * (z[,1] < 1.5) * (z[,2] >= 0.5) * (-5 + x[,1:(p/2)] %*% betas[1:(p/2)] * 2) +
    1 * (z[,1] >= 1.5) * (z[,3] < 0.5) * (5 + x[,1:(p/2)] %*% betas[1:(p/2)] * 0.5) +
    1 * (z[,1] >= 1.5) * (z[,3] >= 0.5) * (10 + x[,1:(p/2)] %*% betas[1:(p/2)] * 0.125) +
    x[,(p/2 + 1):p] %*% betas[(p/2 + 1):p] 
}

# Parameters

load("correlationmatrix.Rdata")  # Should load CorrX
Nsim <- 500
N_vals <- c(250,500)
p_vals <- c(4, 8, 16, 24, 32, 40)
results <- data.frame()
set.seed(4)
betas <- matrix(rlaplace(48, 0, 0.6), ncol = 1)
betas
for (N in N_vals) {
  for (p in p_vals) {
    cat("Running simulations for N =", N, ", p =", p, "\n")
    
    ids <- sample(1:ncol(CorrX), p)
    CorrX1 <- CorrX[ids, ids]
    
    for (sim in 1:Nsim) {
      cat("Simulation ", sim)
      set.seed(sim)
      
      # Simulate Z
      Z <- data.frame(matrix(runif(N * 5), nrow = N))
      Z[,1] <- sample(0:3, N, replace = TRUE)
      Z[,2] <- sample(0:1, N, replace = TRUE)
      colnames(Z) <- paste0("z", 1:5)
      
      # Simulate X and betas
      X <- mvtnorm::rmvnorm(N, sigma = CorrX1)
      colnames(X) <- paste0("x", 1:p)
      
      betas1 <- betas[1:p]
      
      Y <- g(Z, X, betas1) + rnorm(N, 0, 1)
      Y <- as.numeric(Y)
      
      data_train <- cbind(data.frame(Y = Y), Z, X)
      
      # Generate test data
      N_test <- 5000
      Z_test <- data.frame(matrix(runif(N_test * 5), nrow = N_test))
      Z_test[,1] <- sample(0:3, N_test, replace = TRUE)
      Z_test[,2] <- sample(0:1, N_test, replace = TRUE)
      colnames(Z_test) <- paste0("z", 1:5)
      X_test <- mvtnorm::rmvnorm(N_test, sigma = CorrX1)
      colnames(X_test) <- paste0("x", 1:p)
      Y_test <- g(Z_test, X_test, betas1) + rnorm(N_test, 0, 1)
      Y_test <- as.numeric(Y_test)
      
      ## 1. MOB
      mob_formula <- as.formula(paste(
        "Y ~", paste(colnames(X), collapse = " +"), "|", paste(colnames(Z), collapse = " +")
      ))
      mob_fit <- tryCatch({
        lmtree(mob_formula, data = data_train)
      }, error = function(e) NULL)
      
      if (!is.null(mob_fit)) {
        Y_pred_mob <- predict(mob_fit, newdata = cbind(Z_test, X_test))
        r2_mob <- 1 - sum((Y_test - Y_pred_mob)^2) / sum((Y_test - mean(Y_test))^2)
        results <- rbind(results, data.frame(N = N, p = p, Model = "MOB", R2 = r2_mob))
      }
      
      ## 2. fusedTree
      dat <- cbind.data.frame(Y = Y, Z)
      rp <- rpart(Y ~ ., data = dat, control = rpart.control(xval = 5, maxdepth = 4), model = TRUE)
      cp <- rp$cptable[,1][which.min(rp$cptable[,4])]
      Treefit <- prune(rp, cp = cp)
      r2_fused <- NA  # default NA value
      
      tryCatch({
        if (nrow(rp$frame) > 1) {
          foldsTree <- CVfoldsTree(Y = Y, Tree = Treefit, Z = Z, model = "linear", nrepeat = 3, kfold = 5)
          
          optPenalties <- PenOpt(Treefit, X, Y, Z, model = "linear",
                                 lambdaInit = 10, alphaInit = 10,
                                 loss = "loglik", LinVars = FALSE,
                                 folds = foldsTree, multistart = FALSE)
          
          fit <- fusedTree(Treefit, X, Y, Z, LinVars = FALSE, model = "linear",
                           lambda = optPenalties[1], alpha = optPenalties[2])
          
          pred_fused <- predict(fit, newX = X_test, newY = Y_test, newZ = Z_test)$Ypred
          r2_fused <- 1 - sum((Y_test - pred_fused)^2) / sum((Y_test - mean(Y_test))^2)
        }
      }, error = function(e) {
        message("fusedTree model failed: ", e$message)
        mse_fused <<- NA
      })
      # Store result for fusedTree
      results <- rbind(results, data.frame(N = N, p = p, Model = "fusedTree", R2 = r2_fused))
    }
  }
}
save(results, file = "MOB_vs_fusedTree.Rdata")
load("MOB_vs_fusedTree.Rdata")
# Summarize and plot
library(dplyr)
results_avg <- results %>%
  group_by(N, p, Model) %>%
  summarise(mean_MSE = mean(R2, na.rm = TRUE), .groups = "drop")
library(ggplot2); library(ggsci)
results_avg$Model <- ifelse(results_avg$Model == "fusedTree", "FusedTree", results_avg$Model)
name <- "MOB_vs_FusedTree.pdf"
pdf(name, width=7,height=3)
ggplot(results_avg, aes(x = p, y = mean_MSE, color = Model)) +
  geom_line(linewidth= 1.2) +
  geom_point() + scale_x_continuous(breaks = p_vals) +
  facet_wrap(~N, scales = "free_y", labeller = as_labeller(function(x) paste0("N = ", x))) +
  labs(title = "", 
       x = "Number of variables in the leaves (p)", 
       y = expression(R^2)) +
  theme_light() + scale_color_jama()
dev.off()
