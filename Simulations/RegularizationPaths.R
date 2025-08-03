setwd("~/PhD files/R files/ClinOmics/FusedTree_Paper/Simulations")
source('SideFunctions_Update.R')
source('RidgeFunctions_Update.R')
library(rpart); 
library(rpart.plot)
library(Matrix)
library(VGAM)
library(treeClust)
library(ModTools)
library(reshape2); library(ggplot2); library(viridis)
p <- 10
N <- 500
p_Clin <- 5

set.seed(232)
Z <- matrix(runif(N*p_Clin,0,1),nrow = N, ncol = p_Clin) 
colnames(Z) <- paste0("z",seq(1, p_Clin))
Z <- data.frame(Z)
set.seed(232)
X <- mvtnorm::rmvnorm(N, sigma = diag(p))
colnames(X) <- paste0("x", seq(1, p))

###################### Experiment 1 (EFFECT MODIFICATION) ######################
set.seed(444444)
betas <- matrix(rnorm(p,0,5/p),ncol = 1)
betas
Varids <- c(2,6)
betas[Varids]

g <- function (z,x,betas){
  1*(z[,1] < 0.5)*(z[,2] < 0.5)*(-10 + x[,1:p/5] %*% betas[1:p/5]*6) + 
    1*(z[,1] < 0.5)*(z[,2] >= 0.5)*(-5 + x[,1:p/5]%*% betas[1:p/5]*3) +
    1*(z[,1] >= 0.5)*(z[,4] < 0.5)*(5 + x[,1:p/5]%*% betas[1:p/5]*0.5) +
    1*(z[,1] >= 0.5)*(z[,4] >= 0.5)*(10 + x[,1:p/5]%*% betas[1:p/5]*0.2) +
    x[,(p/5+1):p]%*% betas[(p/5+1):p] 
}

g1 <- function (z){
  1*(z[,1] < 0.5)*(z[,2] < 0.5)*(-10) + 
    1*(z[,1] < 0.5)*(z[,2] >= 0.5)*(-5) +
    1*(z[,1] >= 0.5)*(z[,4] < 0.5)*(5) +
    1*(z[,1] >= 0.5)*(z[,4] >= 0.5)*(10)
}

set.seed(23)
Y1 <- g1(z=Z)
Y <- g(z=Z,x=X,betas = betas)+rnorm(N,0,1)
Y <- Y[,1]
Var(Y)


# fitting tree
dat <- cbind.data.frame(Y,Z)
rp <- rpart(Y~.,data = dat, control = rpart.control(xval = 5, minbucket = 20),
            model = T)
cp <- rp$cptable[,1][which.min(rp$cptable[,4])]
Treefit <- prune(rp, cp =cp)
rpart.plot(Treefit, type = 0, extra = 1,nn=F)

NumNod <- length(unique(Treefit$where))
name <- paste("Tree_RegPath","EffMod.pdf", sep = "_")
pdf(name, width = 5, height = 3)
rpart.plot(Treefit, type = 0,extra = 1,nn = T)
dev.off()



Dat <- Dat_Tree(tree = Treefit, X = X, Z = Z, Y = Y, model = "linear", LinVars = F)
X1 <- Dat$Omics; U1 <- Dat$Clinical; Y1 <- Dat$Response

Dg <- .PenMatr(NumNodes = NumNod, p = p)

beg <- NumNod * (Varids) + 1; end <- NumNod * (Varids) + 1 + 3
ids1 <- c(seq(beg[1], end[1]), seq(beg[2], end[2]))


alpha <- seq(0.001, 5000, 1)

lambda <- c(0.0, 1, 10, 100, 500, 5000)
results <- c()
for (i in 1:length(lambda)) {
  Lam <- lambda[i]+0.01
  Betas <- lapply(alpha, function (x) ridgeGLM2(Y = Y1, U = U1, X = as.matrix(X1), 
                                               lambda = Lam, lambdaG = x, Dg = as.matrix(Dg), 
                                               model = "linear"))
  Betas <- do.call(rbind, Betas)
  colnames(Betas) <- c(colnames(U1), colnames(X1))
  Betas <- Betas[, c(3, ids1)]
  
  
  substring(colnames(Betas), 1, 1)[-1] <- "\u03b2"
  colnames(Betas)[1] <- "c6"
  Res <- cbind.data.frame("Alpha" = log(alpha), Betas)
  melt <- melt(Res, id = "Alpha")
  melt <- cbind.data.frame(melt,"Lambda" = rep(paste("Lambda = ",Lam-0.01), nrow(melt)))
  results <- rbind.data.frame(results, melt)
}

unique(results$Lambda)
results$Lambda_f <- factor(results$Lambda)

results$Lambda_f <- factor(results$Lambda, levels = c("Lambda =  0","Lambda =  1","Lambda =  10",
                                                   "Lambda =  100","Lambda =  500","Lambda =  5000"), 
                           ordered = T)
results$Lambda_f <- gsub("Lambda", "\u03bb", results$Lambda_f) 
names(effect_labels)
effect_labels <- c(
  "Effect_clin" = expression(c[6]),
  "Effect_2,4" = expression(beta[2(N4)]),
  "Effect_2,5" = expression(beta[2(N5)]),
  "Effect_2,6" = expression(beta[2(N6)]),
  "Effect_2,7" = expression(beta[2(N7)]),
  "Effect_6,4" = expression(beta[6(N4)]),
  "Effect_6,5" = expression(beta[6(N5)]),
  "Effect_6,6" = expression(beta[6(N6)]),
  "Effect_6,7" = expression(beta[6(N7)])
)
levels(results$variable) <- names(effect_labels)
name <- paste("RegPath", "EffMod", sep = "_")
library(Cairo)
CairoPDF(name, width = 7, height = 7.5)
plot <- ggplot(results, aes(x = Alpha, y = value, group = variable, colour = variable, linetype = variable)) + 
  geom_line() + 
  theme_light() +
  scale_color_viridis(
    discrete = TRUE,
    option = "D",
    labels = effect_labels
  ) +
  scale_linetype_manual(
    values = c("dashed", rep("solid", 4), rep("dotted", 4)),
    labels = effect_labels
  ) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  labs(
    x = expression(log(alpha)),
    y = "Effect estimates"
  )

# Add facet
plot + facet_wrap(. ~ Lambda_f, nrow = 3)
dev.off()
getwd()
