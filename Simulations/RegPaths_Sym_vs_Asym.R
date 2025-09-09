setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
library(fusedTree); library(rpart); library(rpart.plot)
library(tidyverse)
set.seed(444444)
betas <- matrix(rnorm(p,0,2/p),ncol = 1)
betas
g <- function (z,x,betas){
  1*(z[,1] < 0.5)*(z[,2] < 0.5)*(-10 + x[,1:p/5] %*% betas[1:p/5]*6) +
    1*(z[,1] < 0.5)*(z[,2] >= 0.5)*(-5 + x[,1:p/5]%*% betas[1:p/5]*3) +
    1*(z[,1] >= 0.5)*(z[,4] < 0.5)*(5 + x[,1:p/5]%*% betas[1:p/5]*0.5) +
    1*(z[,1] >= 0.5)*(z[,4] >= 0.5)*(10 + x[,1:p/5]%*% betas[1:p/5]*(-3)) +
    x[,(p/5+1):p]%*% betas[(p/5+1):p]
}

p <- 10
N <- 500
p_Clin <- 4

set.seed(232)
Z <- matrix(runif(N*p_Clin,0,1), nrow = N, ncol = p_Clin)
colnames(Z) <- paste0("z", seq(1, p_Clin))
Z <- data.frame(Z)
set.seed(232)
X <- mvtnorm::rmvnorm(N, sigma = diag(p))
colnames(X) <- paste0("x", seq(1, p))
Y <- g(z = Z, x = X, betas = betas) + rnorm(N, 0, 1)
Y <- Y[,1]


dat <- cbind.data.frame(Y,Z)
rp <- rpart(Y~.,data = dat, control = rpart.control(xval = 5, minbucket = 20),
            model = T)
cp <- rp$cptable[,1][which.min(rp$cptable[,4])]
Treefit <- prune(rp, cp =cp)
rpart.plot(Treefit, type = 0, extra = 1,nn=F)

alpha = seq(0.1, 10000, 5)
lambda = c(0.1, 10, 100,1000)

results_list <- list()

for (i in seq_along(lambda)) {
  Lam <- lambda[i]

  # Symmetric fusion
  Betas_sym = lapply(alpha, function(x) {
    fusedTree(Tree = Treefit, Y = Y, X = X, Z = Z,
              lambda = Lam, alpha = x, symFusion = TRUE,
              model = "linear", LinVars = FALSE)$Effects[c(3,5:8)]
  })

  # Asymmetric fusion
  Betas_asym = lapply(alpha, function(x) {
    fusedTree(Tree = Treefit, Y = Y, X = X, Z = Z,
              lambda = Lam, alpha = x, symFusion = FALSE,
              model = "linear", LinVars = FALSE)$Effects[c(3,5:8)]
  })

  # Convert to data frames
  df_sym <- do.call(rbind, Betas_sym) %>%
    as.data.frame() %>%
    mutate(alpha = alpha,
           lambda = Lam,
           symFusion = "Symmetric") %>%
    pivot_longer(cols = everything()[1:5],  # V1 to V4
                 names_to = "effect_index",
                 values_to = "effect_value")

  df_asym <- do.call(rbind, Betas_asym) %>%
    as.data.frame() %>%
    mutate(alpha = alpha,
           lambda = Lam,
           symFusion = "Asymmetric") %>%
    pivot_longer(cols = everything()[1:5],
                 names_to = "effect_index",
                 values_to = "effect_value")

  results_list[[length(results_list) + 1]] <- df_sym
  results_list[[length(results_list) + 1]] <- df_asym
}

# Combine and clean
results_df <- bind_rows(results_list)

# Label effect indices as 5–8
results_df$effect_index <- factor(results_df$effect_index,
                                  labels = paste0("Effect_", c(1,5:8)))

results_df$lambda <- factor(results_df$lambda,
                            levels = sort(unique(results_df$lambda)))
results_df$lambda_lab <- paste0("lambda = ", results_df$lambda)
results_df$lambda_lab <- gsub("lambda", "\u03bb", results_df$lambda_lab)

# Plot
effect_labels <- c(
  "Effect_1" = expression(c[1]),
  "Effect_5" = expression(beta[1(N4)]),
  "Effect_6" = expression(beta[1(N5)]),
  "Effect_7" = expression(beta[1(N6)]),
  "Effect_8" = expression(beta[1(N7)])
)
getwd()
library(Cairo); library(ggsci); library(viridis)
name <- paste("RegPathPower","Sym_vs_Asym", sep = "_")
CairoPDF(name, width = 7,height = 4)
ggplot(results_df, aes(x = alpha, y = effect_value, color = effect_index)) +
  geom_line(linewidth = 1.1) + scale_x_log10() +
  facet_wrap(~ lambda_lab, ncol = 2) +
  scale_color_jama(labels = effect_labels) +
  labs(
    title = "",
    x = expression(log(alpha)),
    y = "Effect Estimate",
    color = "Effect Index"
  ) +
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "right")
dev.off()








# Fill upper triangle and copy to lower (symmetry)
for (m in 1:(M - 1)) {
  for (m2 in (m + 1):M) {
    diff <- base::abs(ranks[m] - ranks[m2])
    rank_diff[m, m2] <- diff
    rank_diff[m2, m] <- diff
  }
}

# Set diagonal to 1
base::diag(rank_diff) <- 0

# Convert to weight matrix
w_matrix <- 1 / (1 + rank_diff)

# Normalize
w_matrix <- w_matrix/base::rowSums(w_matrix)
