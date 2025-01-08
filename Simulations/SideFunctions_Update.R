Dat_Tree <- function(tree,X,Z,Y,model, LinVars=T) {
  #convert data to right format for regressions: X <- diag(X_node1,X_node2, .....,X_nodeNodeNum)
  #                                              Y <- (Y_node1,Y_node2, .....,Y_nodeNodeNum)
  
  # tree: fit of regression tree
  # X: omics covariate design matrix (intercept (=column with 1's) included)
  # Z: clinical covariate matrix (used for tree fitting, so no intercept required!)
  # Y: response vector (should be same length as nrow(X))
  # idVars: Clinical variables not included in tree, will be modeled linearly
  
  
  ### loading R packages
  if (!require(rpart, quietly = T)) {stop("rpart not installed")}
  if (!require(treeClust, quietly = T)) {stop("treeClust not installed")}
  if (!require(Matrix, quietly = T)) {stop("Matrix not installed")}
  
  
  
  
  
  if (class(LinVars) != "logical"){stop("LinVars is not logical, specify as either TRUE or FALSE")}
  
  
  if (ncol(X) == 0 || nrow(X) == 0){stop("X matrix has either zero rows or zero columns")}
  if (ncol(Z) == 0 || nrow(Z) == 0){stop("Z dataframe has either zero rows or zero columns")}
  if (length(Y) == 0){stop("Y vector is empty")}
  
  if (class(X)[1]!="matrix"){stop("X should be specified as a matrix (not a data.frame)")}
  if (class(Z)!="data.frame"){stop("Z should be specified as a data.frame")}
  if(!(class(Y) == "numeric" || class(Y) =="Surv")) {stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  
  nodes <-  row.names(tree$frame)[rpart.predict.leaves(tree,Z)]
  nodes <- as.numeric(nodes)
  Nodenames = paste0("N",sort(unique(as.numeric(row.names(tree$frame)[tree$where]))))
  
  p = ncol(X)
  if (is.null(colnames(X))){colnames(X) = paste0("x",seq(1,ncol(X)))}
  namesX = colnames(X)
  if (is.null(colnames(Z))){colnames(Z) = pasteo("z",seq(1,ncol(Z)))}
  namesZ = colnames(Z)
  
  NumNodes <- length(unique(tree$where))
  if (NumNodes < 2){
    print("Tree has single node, return design matrices")
    return(list(Clinical = model.matrix(~., Z), Omics = X, Response = Y))
  } else {
    ClinIntercepts = model.matrix(~0+factor(nodes, 
                                            levels = sort(unique(as.numeric(row.names(tree$frame)[tree$where])))))
  
    X_tot <- t(Matrix::KhatriRao(t(X),t(ClinIntercepts)))
    colnames(X_tot) <- c(sapply(namesX, function (x) paste0(x,paste0("_",Nodenames))))
    colnames(ClinIntercepts) = Nodenames
    #if (model=="surv"){
      #ClinIntercepts<-ClinIntercepts[,-1]
      #colnames(ClinIntercepts) = Nodenames[-1]
    #} else {
    #  colnames(ClinIntercepts) = append("Intercept",Nodenames[-1])
    #}
    
    
    if (LinVars==T){
      idVars <- which(sapply(Z, is.numeric) & apply(Z,2,function (x) length(unique(x)) > 5))
      nameZ = colnames(Z)[idVars]
      Clinical <- as.matrix(cbind(ClinIntercepts,Z[,idVars]))
      colnames(Clinical)[-(1:NumNodes)] <- nameZ
    } else {
      Clinical <- ClinIntercepts
    }
    return(list(Clinical = Clinical,Omics = X_tot, Response = Y))
  }
}

PenOpt <- function(Tree,X,Y,Z, model=NULL, 
                   lambdaInit=10,alphaInit=10,
                   lambdaMin=10^(-5),alphaMin=10^(-5),
                   folds=.CVfoldsTree(Y=Y, Tree = Tree, model = model,kfold = 5,nrepeat = 1),
                   loss="loglik",
                   maxIter = 100,
                   LinVars=F) {
                   
  
  if (class(Tree) != "rpart"){stop("Tree is not rpart object")}
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  if(!(loss == "loglik" || loss == "sos")){stop("Loss should be specified as loglik or sos")}
  if (ncol(X) == 0 || nrow(X) == 0){stop("Omics covariates not specified")}
  if (ncol(Z) == 0 || nrow(Z) == 0){stop("Clinical covariates  not specified")}
  if (length(Y) == 0){stop("Y vector is empty")}
  
  if(!(class(Y) == "numeric" || class(Y) =="Surv")) {stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  if (class(X)[1]!="matrix"){stop("X should be specified as matrix object (so not a data.frame)")}
  if(class(Z)[1] == "matrix" ){Z <- data.frame(Z)}
  if (class(Z)!="data.frame"){stop("Z should be specified as a data.frame or matrix")}
  
  if (class(LinVars) != "logical"){stop("LinVars is not logical, specify as either TRUE or FALSE")}
  if(class(folds) != "list"){stop("Test fold ids should be collected in list")}
  if(length(folds)<3){stop("Number of folds should be at least 3")}
  if(!all(unlist(lapply(folds,class)) == "integer" | unlist(lapply(folds,class)) == "numeric")){stop("Test fold ids should be specified as integer")}
  if(!all(unlist(folds) <= length(Y))){stop("fold id out of range")}
  
  ## obtaining tree stats
  
  nodes <-  row.names(Tree$frame)[rpart.predict.leaves(Tree,Z)]
  nodes <- as.numeric(nodes)
  names = paste0("N",sort(unique(nodes)))
  
  NumNodes <- length(unique(nodes))
  
  if (NumNodes < 2){
    print("Only single node present in the tree, fitting ridge regression with X penalized and Z unpenalized instead")
    print("Estimating standard penalty lambda for ridge regression")
    
    if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
    U = model.matrix(~.,Z)
    Lam <- optPenaltyGLM.kCVauto2(Y = Y, X = X, U = U, 
                                  lambdaInit = lambdaInit, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)),
                                  model = model, folds = folds, loss = loss,
                                  lambdaMin = lambdaMin,lambdaGmin = alphaMin,
                                  maxIter = maxIter)
    names(Lam) <- "lambda"
    return(Lam)
  }
  
  if (NumNodes > 1){
    Dat = Dat_Tree(tree = Tree, X = X, Y = Y, Z = Z,  LinVars = LinVars, model = model)
   
    X1 = Dat$Omics; Y1 = Dat$Response; U1 = Dat$Clinical; 
    remove(Dat)
    
    NumNod = ncol(X1)/ncol(X)
    Delta = .PenMatr(NumNodes = NumNod, p = ncol(X))
    if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
    if(alphaInit <= 0){stop("initial penalty alphaInit should be larger than zero")}
    if (ncol(Delta) != ncol(X1)){stop("number of columns of penalty matrix does not equal number of columns of design matrix X")}
    optPenalties <- optPenaltyGLM.kCVauto2(Y = Y1, X = X1, U = U1, 
                                           lambdaInit = lambdaInit, lambdaGinit = alphaInit, Dg=Delta,
                                           model = model, folds = folds, loss = loss,
                                           lambdaMin = lambdaMin,lambdaGmin = alphaMin,
                                           maxIter = maxIter)
    names(optPenalties)<-c("lambda","alpha")
    return(optPenalties)
  }
}
    



FusTreeFit <- function(Tree,X,Y,Z,LinVars=F, model, 
                       lambda,alpha, maxIter = 100,minSuccDiff=10^(-10), 
                       dat = T) {
  
  if (class(Tree) != "rpart"){stop("Tree is not rpart object")}
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  
  if(length(lambda) != 1){stop("Lambda should be single numeric value")}
  if(lambda <= 0){stop("penalty lambda should be larger than zero")}
  
  
  if (ncol(X) == 0 || nrow(X) == 0){stop("Omics covariates not specified")}
  if (ncol(Z) == 0 || nrow(Z) == 0){stop("Clinical covariates  not specified")}
  if (length(Y) == 0){stop("Y vector is empty")}
  
  if(!(class(Y) == "numeric" || class(Y) =="Surv")) {stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  if (class(X)[1]!="matrix"){stop("X should be specified as matrix object (so not a data.frame)")}
  if (class(Z)!="data.frame"){stop("Z should be specified as a data.frame")}
  
  
  
  if (class(dat) != "logical"){stop("dat is not logical, specify as either TRUE or FALSE")}
  if (class(LinVars) != "logical"){stop("LinVars is not logical, specify as either TRUE or FALSE")}
  
  
  nodes <-  row.names(Tree$frame)[rpart.predict.leaves(Tree,Z)]
  nodes <- as.numeric(nodes)
  names = paste0("N",sort(unique(nodes)))
  
  NumNodes <- length(unique(nodes))
  
  if (NumNodes < 2){
    print("Number of nodes is 1, fitting ridge regression with X penalized and Z unpenalized instead")
    print(paste("Fitting ridge with standard ridge penalty Lambda = ",round(lambda,1),sep = ""))
    U = model.matrix(~.,Z)
    if (is.null(colnames(U))){colnames(U) = pasteo("u",seq(1,ncol(U)))}
    if(lambda <= 0){stop("penalty lambda should be larger than zero")}
    Fit <- ridgeGLM2(Y = Y, U = U, X = X, 
                     lambda = lambda, lambdaG = 0, Dg = matrix(0, ncol=ncol(X), nrow = ncol(X)), 
                     model=model, maxIter = maxIter, minSuccDiff = minSuccDiff)
    names(Fit)<-c(colnames(U),colnames(X))
    Pars = cbind.data.frame("Model"= model, "LinVar" = LinVars,"Lambda" = lambda)
    if (dat==T){
      return(list(Tree = Tree,Effects = Fit,Omics=X,Clinical=Z,Response = Y,Pars = Pars))
    } else {return(list(Tree = Tree, Effects=Fit,Pars = Pars))}
  }
  
  
  if (NumNodes > 1){
    
    Dat = Dat_Tree(tree = Tree, X = X, Y = Y,Z = Z,  LinVars = LinVars, model = model)
    
    X1 = Dat$Omics; Y1 = Dat$Response; U1 = Dat$Clinical
    remove(Dat)
    if(alpha <= 0){stop("penalty alpha should be larger than zero")}
    if(lambda <= 0){stop("penalty lambda should be larger than zero")}
    NumNod = ncol(X1)/ncol(X)
    Delta = .PenMatr(NumNodes = NumNod, p = ncol(X))
    Fit = ridgeGLM2(Y = Y1, U = U1, X = X1, 
                    lambda = lambda, lambdaG = alpha, Dg = Delta, 
                    model=model, maxIter = maxIter, minSuccDiff = minSuccDiff)
    names(Fit)<-c(colnames(U1),colnames(X1))
    Pars = cbind.data.frame("Model"= model, "LinVar" = LinVars,"Alpha" = alpha,"Lambda" = lambda)
    if (dat==T){
      return(list(Tree= Tree, Effects = Fit,Omics=X1,Clinical=U1,Response = Y1, Pars = Pars))
    } else {return(list(Tree= Tree, Effects = Fit,Pars = Pars))}
  }
}
    


Predictions <- function(fit,newX,newZ,newY, Dat = T,model, Linvars){
  #if (class(Tree) != "rpart"){stop("Tree is not rpart object")}
  
  if (ncol(newX) == 0 | nrow(newX) == 0){stop("Omics matrix has either zero columns or zero rows")}
  if (ncol(newZ) == 0| nrow(newZ) == 0){stop("Clinical matix  has either zero columns or zero rows")}
  if (length(newY) == 0){stop("Y vector is empty")}
  
  if(!(class(newY) == "numeric" || class(newY) =="Surv")) {stop("Y is not a numeric or Surv object. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  if (class(newX)[1]!="matrix"){stop("X should be specified as matrix object (so not a data.frame)")}
  if (class(newZ)!="data.frame"){stop("Z should be specified as a data.frame")}
  if (length(fit)==0){stop("fit object is empty")}
  
  Tree = fit$Tree
  Ests = fit$Effects
  Pars = fit$Pars
  Linvars = Pars$LinVar
  model = Pars$Model
  
  
  Dat_tree_test = Dat_Tree(tree=Tree, X=newX, Y=newY, Z=newZ, LinVars = Linvars, model = model)
  
  
  
  Xte_tree = Dat_tree_test$Omics; Ute_tree = Dat_tree_test$Clinical; Yte_tree=Dat_tree_test$Response
  # check if some covariates need to be removed from omics test data set (if global test says so)
  # check if there are empty nodes for test data. If so remove node specific estimates (omics and intercept) from fit
  ids1 = which(colnames(Xte_tree) %in% names(Ests))
  Xte_tree <- Xte_tree[,ids1]
  
  ids = which(names(Ests) %in% c(colnames(Ute_tree), colnames(Xte_tree)))
  Ests = Ests[ids]
  
  LP = as.numeric(cbind(Ute_tree,Xte_tree) %*% Ests)
  
  if (model=="linear"){
    Ypred=LP
    if (Dat==T){
      return(list(Preds = cbind.data.frame(Resp=Yte_tree,Ypred=Ypred),Fit=Ests,Xtest=Xte_tree,Utest=Ute_tree))
    } else {
      return(list(Preds = cbind.data.frame(Resp=Yte_tree,Ypred=Ypred),Fit=Ests))
    }
  }
  
  if (model=="logistic"){
    Ypred=exp(LP)/(1+exp(LP))
    if (Dat==T){
      return(list(Preds = cbind.data.frame(Resp=Yte_tree,Ypred=Ypred,LinPred = LP),Fit=Ests,Xtest=Xte_tree,Utest=Ute_tree))
    } else {
      return(list(Preds = cbind.data.frame(Resp=Yte_tree,Ypred=Ypred,LinPred = LP),Fit=Ests))
    }
  }
  
  if (model=="surv"){
    Breslow = .breslow(Y=Yte_tree,lp=LP)
    if (Dat==T){
      return(list(Preds = cbind.data.frame(Resp=Yte_tree,BaseCumHaz=Breslow$Ht,BaseHaz=Breslow$ht,LinPred = LP),Fit=Ests,Xtest=Xte_tree,Utest=Ute_tree))
    } else {
      return(list(Preds = cbind.data.frame(Resp=Yte_tree,BaseCumHaz=Breslow$Ht,BaseHaz=Breslow$ht,LinPred = LP),Fit=Ests))
    }
  }
}

Backward_FusedTree <- function(Tree,X,Y,Z,LinVars=T, model, 
                               folds = .CVfoldsTree(Y=Y,Tree=Tree, Z=Z, model=model,kfold=5, nrepeat=1),
                               lambdaInit = 10, alphaInit = 10,lambdaMin = 10^(-5), 
                               alphamin = 10^(-5),
                               loss = "loglik"){
  
  # Test if everything is in right format
  if (class(Tree) != "rpart"){stop("Tree is not rpart object")}
  rpart.plot(Tree, type = 0,extra = 1,nn=T)
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  if (ncol(X) == 0 || nrow(X) == 0){stop("Omics covariates not specified")}
  if (ncol(Z) == 0 || nrow(Z) == 0){stop("Clinical covariates  not specified")}
  if (length(Y) == 0){stop("Y vector is empty")}
  if (class(Y) != "numeric") {stop("Y is not a numeric. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  if (class(X)[1]!="matrix"){stop("X should be specified as matrix object (so not a data.frame)")}
  if (class(Z)!="data.frame"){stop("Z should be specified as a data.frame")}
  if (class(LinVars) != "logical"){stop("LinVars is not logical, specify as either TRUE or FALSE")}
  
  # load required R packages
  if (!require(globaltest, quietly = T)) {stop("globaltest not installed")}
  if (!require(rpart, quietly = T)) {stop("rpart not installed")}
  if (!require(rpart.plot, quietly = T)) {stop("rpart.plot not installed")}
  if (!require(treeClust, quietly = T)) {stop("treeClust not installed")}
  if (!require(Matrix, quietly = T)) {stop("Matrix not installed")}
  
  
  rpart.plot(Tree, type = 0,extra = 1,nn=T)
  nodes <-  row.names(Tree$frame)[rpart.predict.leaves(Tree,Z)]
  nodes <- as.numeric(nodes)
  Nodenames = paste0("N",sort(unique(nodes)))
  p = ncol(X)
  NumNodes <- length(unique(nodes))
  
  if(NumNodes < 2){
    return("Tree has only 1 node, removing omics effects is therefore not required. Please consider function FusTreeFit instead")
  }
  pvals = sapply(sort(unique(nodes)), function(x) globaltest::p.value(globaltest::gt(
    Y[nodes==x],null = rep(1,length(Y[nodes==x])), 
    alternative = data.frame(matrix(X[nodes == x, ],nrow =length(which(nodes==x)),ncol=ncol(X))),
    model = model)))
  names(pvals) <- Nodenames
  
  print("Estimated p values in the nodes are:")
  print(pvals)
  print("Remove omics effects in nodes in order:")
  pvals<-sort(pvals,decreasing = T)
  print(paste(names(pvals),collapse = " -> "))
  
  res <- matrix(NA,ncol=length(pvals)+1,nrow=length(folds))
  
  names <- "None"
  for (l in 1:length(pvals)) {
    names <- append(names,paste(names(pvals)[1:l],collapse = ", "))
  }
  colnames(res)<-names
  rownames(res)<- paste0("Fold ",1:length(folds))
  
  for (i in 1:length(folds)){
    print(paste("Fold", i, sep = " "))
    
    ids= folds[[i]]
    Xtr = X[-ids,]; Xte = X[ids,]
    Ztr = Z[-ids,]; Zte = Z[ids,]
    Ytr = Y[-ids]; Yte = Y[ids]
    remove(ids)
    LogLikFold <- c()
    
    print("Fit Full FusedTree model first")
    Delta = .PenMatr(NumNodes = NumNodes, p = p)
    
    if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
    if(alphaInit <= 0){stop("initial penalty alphaInit should be larger than zero")}
    Dat = Dat_Tree(tree = Tree, X = Xtr, Y = Ytr,Z = Ztr,  LinVars = LinVars)
    X1 = Dat$Omics; Y1 = Dat$Response; U1 = Dat$Clinical
    remove(Dat)
    foldsHyp = .CVfoldsTree(Y1,Tree = Tree,Z=Ztr,model=model,kfold=5,nrepeat=1)
    
    optPenalties <- optPenaltyGLM.kCVauto2(Y = Y1, X = X1, U = U1, 
                                           lambdaInit = lambdaInit, lambdaGinit = alphaInit, 
                                           #lambdaMin = lambdaMin, lambdaGmin = alphamin,
                                           folds = foldsHyp,
                                           Dg=Delta,
                                           model = model,
                                           loss=loss)
    names(optPenalties)<-c("lambda","alpha")
    Fit = ridgeGLM2(Y = Y1, U = U1, X = X1, 
                    lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                    model=model)
    names(Fit)<-c(colnames(U1),colnames(X1))
    Dat_tree_test = Dat_Tree(tree=Tree, X=Xte, Y=Yte, Z=Zte, LinVars = LinVars)
    Xte_tree = Dat_tree_test$Omics; Ute_tree = Dat_tree_test$Clinical; Yte_tree=Dat_tree_test$Response
    remove(Dat_tree_test)
    # check if some covariates need to be removed from omics test data set (if global test says so)
    # check if there are empty nodes for test data. If so remove node specific estimates (omics and intercept) from fit
    ids = which(names(Fit) %in% c(colnames(Ute_tree), colnames(Xte_tree)))
    Ests = Fit[ids]
    
    LP = as.numeric(cbind(Ute_tree,Xte_tree) %*% Ests)
    if (model=="linear"){
      LogLikFold[1] <- sum((Yte_tree-LP)^2)
    }
    if (model=="logistic"){
      LogLikFold[1] <- -.loglikBLMlp(Yte_tree,lp=LP)
    }
    remove(Fit,Ests,LP)  
    
    for (j in 1:length(pvals)) {
      
      EmpNodes = names(pvals)[1:j]
      #names(Fits)[i] <- paste(names(pvals)[1:i],collapse = ",")
      print(paste("Fit FusedTree without omics effects in", paste(names(pvals)[1:j],collapse = ", ")))
      ids = lapply(EmpNodes, function (x) grepl(x,colnames(X1)))
      ids1 = which(Reduce("|",ids))
      X2 <- X1[,-ids1]
      NodNum <- ncol(X2)/p
      
      if (NodNum == 0){
        
        if (model=="linear"){fam = "gaussian"};if (model=="logistic"){fam="binomial"}
        
        
        Fit <- glm(Y1 ~ 0 + ., data = data.frame(U1),family = fam)$coefficients
        Dat_tree_test = Dat_Tree(tree=Tree, X=Xte, Y=Yte, Z=Zte, LinVars = LinVars)
        Ute_tree = Dat_tree_test$Clinical; Yte_tree=Dat_tree_test$Response
        
        ids = which(names(Fit) %in% colnames(Ute_tree))
        Ests = Fit[ids]
        
        LP <- as.numeric(Ute_tree %*% Ests)
        if (model=="linear"){
          LogLikFold[1+j] <- sum((Yte_tree-LP)^2)
        }
        if (model=="logistic"){
          LogLikFold[1+j] <- -.loglikBLMlp(Yte_tree,lp=LP)
        }
        remove(Ests,LP,Fit,ids)
      }
      if (NodNum == 1){
        
        Lam <- optPenaltyGLM.kCVauto2(Y = Y1, X = X2, U = U1, 
                                      lambdaInit = lambdaInit, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X2), nrow = ncol(X2)),
                                      folds = foldsHyp,
                                      model = model, loss=loss)
        Fit <-ridgeGLM2(Y = Y1, U = U1, X = X2, 
                        lambda = Lam, lambdaG = 0, Dg = matrix(0, ncol=ncol(X2), nrow = ncol(X2)), 
                        model=model)
        names(Fit)<-c(colnames(U1),colnames(X2))
        Dat_tree_test = Dat_Tree(tree=Tree, X=Xte, Y=Yte, Z=Zte, LinVars = LinVars)
        Xte_tree = Dat_tree_test$Omics; Ute_tree = Dat_tree_test$Clinical; Yte_tree=Dat_tree_test$Response
        ids1 = which(colnames(Xte_tree) %in% names(Fit))
        Xte_tree <- Xte_tree[,ids1]
        
        ids = which(names(Fit) %in% c(colnames(Ute_tree), colnames(Xte_tree)))
        Ests = Fit[ids]
        LP = as.numeric(cbind(Ute_tree,Xte_tree) %*% Ests)
        
        if (model=="linear"){
          LogLikFold[j+1] <- sum((Yte_tree-LP)^2)
        }
        if (model=="logistic"){
          LogLikFold[j+1] <- -.loglikBLMlp(Yte_tree,lp=LP)
        }
        remove(Ests,LP,Fit,ids,ids1)
        
      }
      if (NodNum > 1){
        
        Delta = .PenMatr(NumNodes = NodNum, p = p)
        if (ncol(Delta) != ncol(X2)){stop("number of columns of penalty matrix does not equal number of columns of design matrix X")}
        if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
        if(alphaInit <= 0){stop("initial penalty alphaInit should be larger than zero")}
        optPenalties <- optPenaltyGLM.kCVauto2(Y = Y1, X = X2, U = U1, 
                                               lambdaInit = lambdaInit, lambdaGinit = alphaInit,
                                               folds = foldsHyp,
                                               Dg=Delta,model = model,loss=loss)
        
        
        Fit = ridgeGLM2(Y = Y1, U = U1, X = X2, 
                        lambda = optPenalties[1], lambdaG = optPenalties[2], Dg = Delta, 
                        model=model)
        names(Fit)<-c(colnames(U1),colnames(X2))
        Dat_tree_test = Dat_Tree(tree=Tree, X=Xte, Y=Yte, Z=Zte, LinVars = LinVars)
        Xte_tree = Dat_tree_test$Omics; Ute_tree = Dat_tree_test$Clinical; Yte_tree=Dat_tree_test$Response
        ids1 = which(colnames(Xte_tree) %in% names(Fit))
        Xte_tree <- Xte_tree[,ids1]
        ids = which(names(Fit) %in% c(colnames(Ute_tree), colnames(Xte_tree)))
        Ests = Fit[ids]
        LP = as.numeric(cbind(Ute_tree,Xte_tree) %*% Ests)
        
        if (model=="linear"){
          LogLikFold[j+1] <- sum((Yte_tree-LP)^2)
        }
        if (model=="logistic"){
          LogLikFold[j+1] <- -.loglikBLMlp(Yte_tree,lp=LP)
        }
        remove(Ests,LP,Fit,ids,ids1)
      }
      remove(X2)
    }
    res[i,]<-LogLikFold
    remove(LogLikFold,Xtr,Xte,Ztr,Zte,Ytr,Yte,ids1,foldsHyp,X1,Y1,U1)
  }
  
  ave = colMeans(res)
  win = which.min(ave)
  print(paste("Winner is with",colnames(res)[win],"removed",sep = " "))
  print("Fit winner on all data")
  Dat = Dat_Tree(tree = Tree, X = X, Y = Y,Z = Z,  LinVars = LinVars)
  X1 = Dat$Omics; Y1 = Dat$Response; U1 = Dat$Clinical
  remove(Dat)
  foldsHyp = .CVfoldsTree(Y1,Tree = Tree,Z=Z,model=model,kfold=5,nrepeat=1)
  
  if(win>1){
    EmpNodes = names(pvals)[1:(win-1)]
    ids = lapply(EmpNodes, function (x) grepl(x,colnames(X1)))
    ids1 = which(Reduce("|",ids))
    X1 <- X1[,-ids1]
    remove(ids,ids1)
  } else{
    X1 <- X1
  }
  NodNum <- ncol(X1)/p
  
  if (NodNum == 0){
    print("Omics covariates in all nodes set to zero")
    print("Only estimate clinical effects by unpenalized regression")
    if (model=="linear"){fam = "gaussian"};if (model=="logistic"){fam="binomial"}
    Fit <- glm(Y1 ~ 0 + ., data = data.frame(U1),family = fam)$coefficients
    Lam = "Not required"
    Alf = "Not required"
  }
  if (NodNum == 1){
    print(paste("Only 1 node with omics effect, fusion penalty alpha is not necessary. Only use standard ridge penalty lambda"))
    optPenalties <- optPenaltyGLM.kCVauto2(Y = Y1, X = X1, U = U1, 
                                           lambdaInit = lambdaInit, lambdaGinit = 0, Dg = matrix(0, ncol=ncol(X1), nrow = ncol(X1)),
                                           folds = foldsHyp,
                                           model = model, loss=loss)
    
    Lam <- optPenalties[1]
    Alf <- "Not required"
    Fit <- ridgeGLM2(Y = Y1, U = U1, X = X1, 
                    lambda = Lam, lambdaG = 0, Dg = matrix(0, ncol=ncol(X1), nrow = ncol(X1)), 
                    model=model)
    names(Fit)<-c(colnames(U1),colnames(X1))
  }
  if (NodNum > 1){
    Delta = .PenMatr(NumNodes = NodNum, p = p)
    if (ncol(Delta) != ncol(X1)){stop("number of columns of penalty matrix does not equal number of columns of design matrix X")}
    if(lambdaInit <= 0){stop("initial penalty lambdaInit should be larger than zero")}
    if(alphaInit <= 0){stop("initial penalty alphaInit should be larger than zero")}
    optPenalties <- optPenaltyGLM.kCVauto2(Y = Y1, X = X1, U = U1, 
                                           lambdaInit = lambdaInit, lambdaGinit = alphaInit,
                                           folds = foldsHyp,
                                           Dg=Delta,model = model,loss=loss)
    Lam <- optPenalties[1]
    Alf <- optPenalties[2]
    Fit <- ridgeGLM2(Y = Y1, U = U1, X = X1, 
                    lambda = Lam, lambdaG = Alf, Dg = Delta, 
                    model=model)
    names(Fit)<-c(colnames(U1),colnames(X1))
    
  }
  Pars = cbind.data.frame("Model"= model, "LinVar" = LinVars,"Lambda" = Lam,"Alpha" = Alf)
  return(list(CVPerformance = res,
              Tree = Tree, 
              Effects = Fit,
              Omics = X, Clinical = Z, Response = Y,
              Pars = Pars))
}





BaggedFusedTree <- function(X,Y,Z,LinVars=T, model, 
                            nBag = 10,
                            lambdaInit = 10, alphaInit = 10,
                            lambdaMin = 10^(-5),alphamin = 10^(-5),
                            loss = "loglik"){
  # Test if everything is in right format
  
  
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  if(!(loss == "loglik" || loss == "sos")){stop("Loss should be specified as loglik or sos")}
  if (ncol(X) == 0 || nrow(X) == 0){stop("Omics covariates not specified")}
  if (ncol(Z) == 0 || nrow(Z) == 0){stop("Clinical covariates  not specified")}
  if (length(Y) == 0){stop("Y vector is empty")}
  if (class(Y) != "numeric") {stop("Y is not a numeric. If Y is binary please specify it as numeric vector coded with 0 and 1")}
  if (class(X)[1]!="matrix"){stop("X should be specified as matrix object (so not a data.frame)")}
  if (class(Z)!="data.frame"){stop("Z should be specified as a data.frame")}
  if (class(LinVars) != "logical"){stop("LinVars is not logical, specify as either TRUE or FALSE")}
  
  # load required R packages
  if (!require(globaltest, quietly = T)) {stop("globaltest not installed")}
  if (!require(rpart, quietly = T)) {stop("rpart not installed")}
  if (!require(rpart.plot, quietly = T)) {stop("rpart.plot not installed")}
  if (!require(treeClust, quietly = T)) {stop("treeClust not installed")}
  if (!require(Matrix, quietly = T)) {stop("Matrix not installed")}
  
  Fits <- list()
  
  for (i in 1:nBag) {
    print(paste("Bootstrap sample ",i))
    
    # specify bootstrapped data
    ids=sample(1:length(Y),length(Y),replace = T)
    Yboot <- Y[ids]; Xboot <- X[ids,]; Zboot <- Z[ids,]
    
    # Fit tree
    dat=cbind.data.frame(Yboot,Zboot)
    Treefit <- rpart(Yboot~.,data = dat, control = rpart.control(xval =10, minbucket = 25,cp=0),
                model = T)
    
    rpart.plot(Treefit, type = 0,extra = 1,nn=T)
    Nodes <-  row.names(Treefit$frame)[rpart.predict.leaves(Treefit,Z)]
    Nodes <- as.numeric(Nodes)
    NumNod <- length(unique(Nodes))
    # Fit FusedTree
    optPenalties <- PenOpt(Tree = Treefit,X = Xboot, Y = Yboot, Z = Zboot, lambdaInit = lambdaInit, alphaInit = alphaInit,
                           model = model, loss = loss,
                           LinVars = LinVars)
    
    Fit = FusTreeFit(Tree = Treefit,X = Xboot,Y = Yboot,Z = Zboot,
                     LinVars = LinVars, model = model, 
                     lambda = optPenalties[1], alpha = optPenalties[2],
                     dat=T)
    Fits[[i]]<-Fit
    names(Fits)[i]<- paste("Bootstrap ",i)
    remove(Fit,optPenalties,Treefit,NumNod,Nodes,Xboot,Yboot,Zboot)
  }
  #aggregate results
  return(Fits) 
  #list with each element the FusedTree fit corresponding to the bootstrapped sample, 
  # bootstrapped data is also found for corresponding Fit
}

.EmptyNodes <- function(Tree,Ztest){
  
  NodesTr <-  row.names(Tree$frame)[Tree$where]
  NodesTr <- as.numeric(NodesTr)
  NodesTe <-  row.names(Treefit$frame)[rpart.predict.leaves(Tree,Ztest)]
  NodesTe <- as.numeric(NodesTe)
  EmpNodes = setdiff(NodesTr,NodesTe)
  if(length(EmpNodes)==0) {return(EmpNodes)}
 
  Empty = paste0("N",EmpNodes,sep = "")
  
  return(Empty)
}


.PenMatr <- function(NumNodes,p){
  if (!require(Matrix, quietly = T)) {stop("Matrix not installed")}
  block=diag(1,nrow = NumNodes,ncol = NumNodes)-1/NumNodes*matrix(1,nrow = NumNodes,ncol = NumNodes)
  Omega=Matrix::kronecker(Diagonal(p,1), block)
  return(Omega)
}

.EigPenmatr <- function(NumNodes,p){
  if (!require(Matrix, quietly = T)) {stop("Matrix not installed")}
  Eig_base = eigen(.PenMatr(NumNodes,1))$vectors
  
  Eigvecs = Matrix::kronecker(Diagonal(p,1),Eig_base)
  diags = rep(c(rep(1,NumNodes-1),0),p)
  return(list(values=diags,vectors=Eigvecs))
}

.SplitVariables <- function(Tree,Z) {
  rules = unlist(Rules(Tree,leafonly = T)$path)
  Vars <- lapply(colnames(Z), function (x) sum(grepl(x,unlist(rules))))
  names(Vars) <- colnames(Z)
  Vars <- unlist(Vars,use.names = T)
  ids = which(Vars==0)
}


.CVfoldsTree <- function(Y,Tree,Z,model=NULL,kfold=5,nrepeat=1){ #response is required for balanced CV
  #response: response vector, length N
  #model: "logistic", "linear", etc
 
  #kfold: scalar, the number of folds
  
  #nrepeat: number of repeats of the CV
  #Output: list object with kfold elements containing the sample indices of the left-out samples per fold
  if (!require(splitTools, quietly = T)) {stop("splitTools not installed")}
  if (!require(treeClust, quietly = T)) {stop("treeClust not installed")}
  
  
  
  if (class(Tree) != "rpart"){stop("Tree is not rpart object")}
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  
  Nodes <-  row.names(Tree$frame)[rpart.predict.leaves(Tree,Z)]
  Nodes <- factor(Nodes)
  if (model=="logistic"){
    StratDat <- cbind.data.frame(Resp=factor(Y), Node = Nodes) # for binary Y stratify on tree cluster and outcome
    Strat <- multi_strata(StratDat, strategy = "interaction")
    folds = create_folds(Strat, k = kfold, invert = T, m_rep = nrepeat)
  }
  if (model=="linear"){
    StratDat <- cbind.data.frame(Node = Nodes) # for continuous Y only stratify on tree cluster
    Strat <- multi_strata(StratDat, strategy = "interaction")
    folds = create_folds(Strat, k = kfold, invert = T, m_rep = nrepeat)
  }
  if (model=="surv"){
    StratDat <- cbind.data.frame(Resp = factor(Y[,2]),Node = Nodes) # for survival balance the number of events
    Strat <- multi_strata(StratDat, strategy = "interaction")
    folds = create_folds(Strat, k = kfold, invert = T, m_rep = nrepeat)
  }
  
  
  return(folds)
}


  
  
  
  



