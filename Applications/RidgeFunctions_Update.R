########## Fitting Ridge Regression ###########
###############################################

ridgeGLM2 <- function(Y, X, U, 
                      lambda, lambdaG=0, 
                      Dg=matrix(0, ncol=ncol(X), nrow=ncol(X)), 
                      model, 
                      minSuccDiff=10^(-10), maxIter=100){
  
  ## ---------------------------------------------------------------------
  ## Ridge estimation of the regression parameter of the 
  ## generalized linear model (GLM)
  ## ---------------------------------------------------------------------
  ## Arguments 
  ## Y           : A numeric being the response vector.
  ## X           : A matrix, the design matrix. The number of rows should
  ##               match the number of elements of Y.
  ## lambda      : A numeric, the ridge penalty parameter.
  
  ## model       : A character indicating which generalized linear model
  ##               model instance is to be fitted.
  ## minSuccDiff : A numeric, the minimum distance between two successive 
  ##               estimates to be achieved.
  ## maxIter     : A numeric specifying the maximum number of iterations.
  ## ---------------------------------------------------------------------
  ## Value:
  ## The ridge estimate of the regression parameter, a numeric.
  ## ---------------------------------------------------------------------
  ## Authors      : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------
  if (ncol(X) == 0 | nrow(X)==0) {stop("Omics matrix has zero columns or rows")}
  if (ncol(U) == 0 | nrow(U)==0) {stop("Clinical matrix has zero columns or rows")}
  if (lambda<=0) {stop("lambda should be positive")}
  if (length(Y) == 0){stop("Response Y is empty")}
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  
  
  if (model == "linear"){
    if (length(unique(Y)) < 3) {stop("Model=linear, but Y has maximally 2 distinct values")}
    return(.ridgeLM(Y, X, U, lambda, lambdaG, Dg))
  }
  
  if (model == "logistic"){
    if(!all(Y==1 | Y==0)) {stop("Logistic model, specify binary response as numeric coded with 0 and 1")}
    return(.ridgeBLM(Y, X, U, lambda, lambdaG, Dg, 
                     minSuccDiff, maxIter))
  }
  
  if (model == "surv"){
    if(class(Y)!="Surv"){stop("Survival Model, please specify Y as survival object")}
    return(.ridgeSurv(Y,X,U,lambda,lambdaG, Dg, minSuccDiff,maxIter))
  }
}


optPenaltyGLM.kCVauto2 <- function(Y, 
                                   X, 
                                   U, 
                                   lambdaInit,
                                   lambdaGinit = 0,
                                   Dg = matrix(0, ncol=ncol(X), nrow=ncol(X)),
                                   model = "linear", 
                                   folds,
                                   loss="loglik",
                                   lambdaMin = 10^(-5), 
                                   lambdaGmin = 10^(-5),
                                   minSuccDiff = 10^(-5), 
                                   maxIter = 100){
  
  ## ---------------------------------------------------------------------
  ## Function finds the optimal penalty parameter of the ridge 
  ## regression estimator of the generalized linear model parameter. The 
  ## optimum is defined as the minimizer of the cross-validated loss 
  ## associated with the estimator. 
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y           : response vector
  ## X           : design matrix
  ## lambdaInit  : initial guess of the optimal penalty parameter
  ## folds       : list with folds
  ## stratified  : boolean, should the data by stratified such that 
  ##               range of the response is comparable among folds
  
  ## model       : A character indicating which generalized linear model
  ##               model instance is to be fitted.
  ## loss        : loss to be used in CV, either "sos" (sum-of-squares) or 
  ##              "loglik" (loglikelihood)
  ## minSuccDiff : A numeric, the minimum distance between two successive 
  ##               estimates to be achieved.
  ## maxIter     : A numeric specifying the maximum number of iterations.
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the optimal cross-validated penalty parameter of the
  ## ridge estimator of the generalized linear regression model parameter
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  ####----------------------------------------------------------------###
  ####----------------------------------------------------------------###
  if (!require(subplex, quietly = T)) {stop("subplex not installed, please install")}
  if (!require(nloptr, quietly = T)) {stop("nloptr not installed, please install")}
  
  if (ncol(X) == 0 | nrow(X)==0) {stop("Omics matrix has zero columns or rows")}
  if (ncol(U) == 0 | nrow(U)==0) {stop("Clinical matrix has zero columns or rows")}
  
  if (length(Y) == 0){stop("Response Y is empty")}
  if(!(model == "linear" || model =="logistic" || model == "surv")){stop("Model should be specified as linear, logistic, or surv")}
  if (lambdaInit<=0) {stop("initial lambda should be positive")}
  
  
  
  if (((max(abs(Dg)) == 0) | lambdaGinit == 0) && (ncol(U) != 0)){
    if (model == "linear"){
      if (ncol(X) >= nrow(X)){ 
        lambdasOpt <- optim(lambdaInit, 
                            .kcvLMloss_PgeqN_UP_noGP,
                            Y=Y,
                            X=X,
                            U=U,
                            folds=folds,
                            method="Brent",
                            lower=lambdaMin,
                            upper=10^(10),
                            loss=loss)$par
        return(lambdasOpt)
      }
      if (ncol(X) < nrow(X)){ 
        lambdasOpt <- optim(lambdaInit, 
                            .kcvLMloss_PleqN_UP_GPandNoGP,
                            Y=Y,
                            X=X,
                            U=U,
                            Dg=lambdaGinit*Dg,
                            folds=folds,
                            method="Brent",
                            lower=lambdaMin,
                            upper=10^(10),
                            loss=loss)$par
        return(lambdasOpt)
      }
    }
    if (model == "logistic"){
      if(!all(Y==1 | Y==0)) {stop("Logistic model, specify binary response as numeric coded with 0 and 1")}
      lambdasOpt <- optim(log(lambdaInit),
                          .kcvBLMloss_UP_noGP,
                          Y=Y, 
                          X=X, 
                          U=U,
                          folds=folds,
                          maxIter=maxIter,
                          minSuccDiff=minSuccDiff,
                          method="Brent",
                          lower = log(0.1), upper = log(10^10),
                          control = list(reltol=1e-5, trace = 1, parscale =log(2))
                          )$par
      return(exp(lambdasOpt))
    }
    if (model == "surv"){
      if(class(Y) != "Surv"){stop("Response is not Surv object")}
      if(loss != "loglik"){stop("For survival model, only loglik loss is implemented")}
      
      lambdasOpt <- optim(log(lambdaInit),
                          .kcvSurvloss_UP_noGP,
                          Y=Y, 
                          X=X, 
                          U=U,
                          folds=folds, 
                          maxIter=maxIter,
                          minSuccDiff=minSuccDiff,
                          method = "Brent",
                          lower = log(0.1), upper = log(10^10),
                          control = list(reltol=1e-5, trace = 1, parscale =log(2))
                          )$par
      return(exp(lambdasOpt))
    }
  }
  
  if (((max(abs(Dg)) != 0) && lambdaGinit != 0) && (ncol(U) != 0)){
    if (model == "linear"){
      if (length(unique(Y)) < 3) {stop("Linear model, but Y has maximally 2 distinct values")}
      if (ncol(X) >= nrow(X)){ 
        NumNodes = length(which(Dg[1,]!=0))
        Dg <- .EigPenmatr(NumNodes = NumNodes, p = ncol(X)/NumNodes) #fast eigendecomposition
        lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
                                  .kcvLMloss_PgeqN_UP_GP,
                                  grad = NULL,
                                  Y=Y,
                                  X=X %*% Dg$vectors,
                                  U=U,
                                  Ds=Dg$values,
                                  folds=folds,
                                  ui = diag(2),
                                  ci=c(lambdaMin, lambdaGmin),
                                  loss=loss,
                                  outer.eps = 10e-5,
                                  control = list(reltol=10e-5))$par
        return(lambdasOpt)
      }
      if (ncol(X) < nrow(X)){ 
        lambdasOpt <- constrOptim(c(lambdaInit, lambdaGinit),
                                  .kcvLMloss_PleqN_UP_GPandNoGP,
                                  grad = NULL,
                                  Y=Y,
                                  X=X,
                                  U=U,
                                  Dg=Dg,
                                  folds=folds,
                                  ui=diag(2),
                                  ci=c(lambdaMin, lambdaGmin),
                                  loss=loss,
                                  outer.eps = 10e-5,
                                  control = list(reltol=10e-5))$par
        return(lambdasOpt)
      }
    }
    if (model == "logistic"){
      if(class(Y) == "Surv"){stop("Response is survival object, while logistic model is specified")}
      if(!all(Y==1 | Y==0)) {stop("Logistic model, specify binary response as numeric coded with 0 and 1")}
      if(loss != "loglik"){stop("For logistic model, only loglik loss is implemented")}
      lambdasOpt <- optim(log(c(lambdaInit, lambdaGinit)),
                                .kcvBLMloss_UP_GP,
                                Y=Y,
                                X=X,
                                U=U,
                                Dg=Dg,
                                folds=folds,
                                maxIter=maxIter,
                                minSuccDiff=minSuccDiff,
                                method = "Nelder-Mead",
                                control = list(reltol=1e-5, parscale = c(log(2),log(2)))
                                )$par
      return(exp(lambdasOpt))
    }
    if (model == "surv"){
      if(class(Y) != "Surv"){stop("Response is not Surv object")}
      if(loss != "loglik"){stop("For survival model, only loglik loss is implemented")}
      lambdasOpt <- optim(log(c(lambdaInit, lambdaGinit)),
                          .kcvSurvloss_UP_GP,
                          Y=Y,
                          X=X,
                          U=U,
                          Dg=Dg,
                          folds=folds,
                          maxIter=maxIter,
                          minSuccDiff=minSuccDiff,
                          method = "Nelder-Mead",
                          control = list(reltol=1e-5, parscale = c(log(2),log(2)))
                          )$par
      return(exp(lambdasOpt))
    }
  }
}


.kcvLMloss_PgeqN_UP_noGP <- function(lambda, Y, X, U, folds, loss){
  
  ## ---------------------------------------------------------------------
  ## Internal function yields the cross-validated loglikelihood of the 
  ## ridge regression estimator with a two-dimensional covariate layout. 
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## lambdas : penalty parameter vector
  ## Y       : response vector
  ## X       : design matrix multiplied by the eigenvector matrix of the
  ##           nonnegative definite matrix that specifies the structure 
  ##           of spatial fused ridge penalty.
  ## U       : design matrix of the unpenalized covariates
  
  ## Ds      : nonnegative eigenvalues of the nonnegative definite matrix 
  ##           that specifies the structure of spatial fused ridge penalty.
  ## folds   : list with folds
  ## loss    : character, either 'loglik' of 'sos', specifying the loss
  ##           criterion to be used in the cross-validation
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # evaluate loss per fold
  cvLoss <- 0
  
  # calculation prior to cv-loop
  XXT     <- tcrossprod(X) / lambda
  
  
  for (k in 1:length(folds)){
    # evaluate the estimator of the unpenalized low-dimensional 
    # regression parameter on the training samples
    if (all(dim(U) > 0)){
      gHat <- solve(crossprod(U[-folds[[k]],,drop=FALSE], 
                              solve(XXT[-folds[[k]], -folds[[k]]] + 
                                      diag(nrow(XXT)-length(folds[[k]])), 
                                    U[-folds[[k]],,drop=FALSE])),
                    crossprod(U[-folds[[k]],,drop=FALSE], 
                              solve(XXT[-folds[[k]], -folds[[k]]] + 
                                      diag(nrow(XXT)-length(folds[[k]])), 
                                    Y[-folds[[k]]])))
      gHat <- as.numeric(gHat)
    }
    
    # evaluate the linear predictor on the left-out samples
    if (all(dim(U) > 0)){
      Yhat <- tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
                                 XXT[-folds[[k]], -folds[[k]]], 
                               Y[-folds[[k]]] - 
                                 as.numeric(crossprod(t(U[-folds[[k]],,drop=FALSE]), gHat))), 
                         XXT[folds[[k]], -folds[[k]], drop=FALSE])
    } else {
      Yhat <- tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
                                 XXT[-folds[[k]], -folds[[k]]],
                               Y[-folds[[k]]]), 
                         XXT[folds[[k]], -folds[[k]], drop=FALSE])
    }
    
    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
        (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) / 
               length(folds[[k]])) + 1) / 2
    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
    }
  }
  
  # average over the folds
  return(cvLoss / length(folds))
}

.kcvLMloss_PleqN_UP_GPandNoGP <- function(lambdas, Y, X, U, Dg, folds, loss){
  
  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the linear regression
  ## model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # evaluate loss per fold
  cvLoss <- 0
  for (k in 1:length(folds)){
    # evaluate regression estimate
    bHat <- ridgeGLM(Y[-folds[[k]]], X[-folds[[k]],,drop=FALSE],
                     U[-folds[[k]],,drop=FALSE], Dg=Dg,
                     lambda =lambdas[1], lambdaG=lambdas[2],
                     model="linear")
    
    # evaluate linear predictor		                 	
    Yhat <- X[folds[[k]],,drop=FALSE] %*% bHat[-c(1:ncol(U))] +
      U[folds[[k]],,drop=FALSE] %*% bHat[c(1:ncol(U))]
    
    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
        (log(2 * pi * sum((Y[folds[[k]]] - as.numeric(Yhat))^2) / 
               length(folds[[k]])) + 1) / 2
    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - as.numeric(Yhat))^2)
    }
  }
  
  # average over the folds
  return(cvLoss / length(folds))
}

.kcvBLMloss_UP_noGP <- function(lambda, Y, X, U, folds, minSuccDiff, maxIter){
  
  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------
  
  # initiate
  cvLoss  <- 0
  lambda <- exp(lambda) # positivity constraint
  
  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){ 
    XXT     <- tcrossprod(X)/lambda
  }
  
  
  for (k in 1:length(folds)){
    
    
    # evaluate the initial linear predictor
    lpOld <- rep(log(mean(Y[-folds[[k]]])/(1-mean(Y[-folds[[k]]]))), length(Y[-folds[[k]]]))
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    
    loglikPrev <- .loglikBLMlp(Y=Y[-folds[[k]]], lp = lpOld) #initial loglikelihood (penalty is zero)
    
    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lpOld))
      W0    <- Ypred * (1 - Ypred)
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <- 
          .Machine$double.eps
      }
      if (ncol(X) >= nrow(X)){ 
        # adjusted response
        Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0
        
        # evaluate unpenalized low-dim regression estimator
        diag(XXT)[-folds[[k]]] <- diag(XXT)[-folds[[k]]] + 1/W0
        slh        <- solve(XXT[-folds[[k]], -folds[[k]]], 
                            Z)
        gHat       <- solve(crossprod(U[-folds[[k]],], 
                                      solve(XXT[-folds[[k]], -folds[[k]]], 
                                            U[-folds[[k]],])),
                            crossprod(U[-folds[[k]],], slh))
        Ug         <- tcrossprod(U, t(gHat))
        
        # evaluate linear predictor without evaluating the estimator
        slh        <- solve(XXT[-folds[[k]], -folds[[k]], drop=FALSE], 
                            Z - Ug[-folds[[k]],])
        diag(XXT)[-folds[[k]]] <- diag(XXT)[-folds[[k]]] - 1/W0				
        lpAll      <- crossprod(XXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- sum(lpAll[-folds[[k]]] * slh) / 2
        lpAll      <- lpAll + Ug
      }
      if (ncol(X) < nrow(X)){ 
        # adjusted response
        Z <- W0*lpOld + (Y[-folds[[k]]] - Ypred)
        
        # evaluate subexpressions of the estimator
        XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*")) 
        tUTX  <- crossprod(X[-folds[[k]], , drop=FALSE], 
                           sweep(U[-folds[[k]], , drop=FALSE], 1, W0, FUN="*"))
        XTZ   <- crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <- crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <- crossprod(sweep(U[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
        diag(XTXpD) <- diag(XTXpD) + lambda
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
                                            crossprod(X[-folds[[k]], , drop=FALSE], 
                                                      sweep(U[-folds[[k]], , drop=FALSE], 
                                                            1, W0, FUN="*")))))))
        Ug <- as.numeric(tcrossprod(gHat, U))
        
        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld + 
                      (Y[-folds[[k]]] - Ypred) - W0*Ug[-folds[[k]]]))
        bHat    <- qr.solve(XTXpD, XWZpT, 
                            tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug				
        penalty <- 0.5 * lambda * sum((bHat)^2) 
      }
      
      # split linear predictor by fold
      lpOld      <- lpAll[-folds[[k]]]
      lpNew      <- lpAll[ folds[[k]]]
      
      
      # compute new likelihood
      loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty
      
      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }
      
      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        #print("Singluar Fit")
        lpNew=rep(0,nrow(Y[folds[[k]]]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        #print("Succesful Fit")
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lpAll
      }
    }
    cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
  }
  # average over the folds
  return(-cvLoss / length(folds))
}

.kcvSurvloss_UP_noGP <- function(lambda, Y, X, U, folds, minSuccDiff, maxIter){
  
  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # initiate
  RespTot <- c()
  LPTot <- c()
  lambda <- exp(lambda) # to ensure positive lambda
  
  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){ 
    XXT     <- tcrossprod(X)/lambda
  }
  
  
  for (k in 1:length(folds)){
    
    # initialization
    penalty <- 0
    lpPrev <- rep(0, length(Y[,2]))
    lpOld <- rep(0, length(Y[-folds[[k]],2]))
    loglikPrev <- .loglikSurv(Y=Y[-folds[[k]],], lp=lpOld) - penalty
    Yev = Y[-folds[[k]],2]
    
    
    for (iter in 1:maxIter){
      # calculate the weights
      H0 <- .breslow(Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]),lpOld)[,2]
      Ypred <- as.numeric(H0 * exp(lpOld)) 
      W0 <- Ypred 
      if (min(W0) <= 2*.Machine$double.eps){
        W0[which(W0 < 2*.Machine$double.eps)] <- 
          2*.Machine$double.eps
      }
      
      if (ncol(X) >= nrow(X)){ 
        # adjusted response
        Z <- lpOld + (Yev - Ypred)/W0
        
        # evaluate unpenalized low-dim regression estimator
        diag(XXT)[-folds[[k]]] <- diag(XXT)[-folds[[k]]] + 1/W0
        slh        <- solve(XXT[-folds[[k]], -folds[[k]]], 
                            Z)
        gHat       <- solve(crossprod(U[-folds[[k]],], 
                                      solve(XXT[-folds[[k]], -folds[[k]]], 
                                            U[-folds[[k]],])),
                            crossprod(U[-folds[[k]],], slh))
        Ug         <- tcrossprod(U, t(gHat))
        
        # evaluate linear predictor without evaluating the estimator
        slh        <- solve(XXT[-folds[[k]], -folds[[k]], drop=FALSE], 
                            Z - Ug[-folds[[k]],])
        diag(XXT)[-folds[[k]]] <- diag(XXT)[-folds[[k]]] - 1/W0				
        lpAll      <- crossprod(XXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- sum(lpAll[-folds[[k]]] * slh) / 2
        lpAll      <- lpAll + Ug
      }
      if (ncol(X) < nrow(X)){ 
        # adjusted response
        Z <- W0*lpOld + (Yev - Ypred)
        
        # evaluate subexpressions of the estimator
        XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*")) 
        tUTX  <- crossprod(X[-folds[[k]], , drop=FALSE], 
                           sweep(U[-folds[[k]], , drop=FALSE], 1, W0, FUN="*"))
        XTZ   <- crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <- crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <- crossprod(sweep(U[-folds[[k]], , drop=FALSE], 1, sqrt(W0), FUN="*"))
        diag(XTXpD) <- diag(XTXpD) + lambda
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
                                            crossprod(X[-folds[[k]], , drop=FALSE], 
                                                      sweep(U[-folds[[k]], , drop=FALSE], 
                                                            1, W0, FUN="*")))))))
        Ug <- as.numeric(tcrossprod(gHat, U))
        
        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld + 
                      (Yev - Ypred) - W0*Ug[-folds[[k]]]))
        bHat <- solve(XTXpD,XWZpT)
        #bHat    <- qr.solve(XTXpD, XWZpT, 
        #                    tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug				
        penalty <- 0.5 * lambda * sum((bHat)^2)
      }
      # split linear predictor by fold
      lpOld <- lpAll[-folds[[k]]]
      lpNew <- lpAll[folds[[k]]]
      
      # compute new likelihood
      loglik <- .loglikSurv(Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]), lpOld) - penalty
      
      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }
      
      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        print("Singluar Fit")
        lpNew=rep(0,nrow(Y[folds[[k]],]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        #print("Succesful Fit")
        break
      } else {
        loglikPrev <- loglik
        lpPrev<-lpAll
      }
      if(iter==maxIter){print("Not converged yet, please increase maxIter")}
    }
    
    RespTot <- base::append(RespTot,Y[folds[[k]],])
    LPTot <- append(LPTot,as.numeric(lpNew))
  }
  
  RespTot1 <- Surv(time=RespTot[,1], event = RespTot[,2])
  CVLoss = .loglikSurv(Y=RespTot1,lp=LPTot)
  
  # average over the folds
  return(-CVLoss/length(folds))
}

.kcvSurvloss_UP_noGP1<- function(lambda, Y, X, U, folds, minSuccDiff, maxIter){
  
  RespAll <- c()
  LPAll <- c()
  lambda<-exp(lambda)
  
  
  for (k in 1:length(folds)){
    ids=folds[[k]]
    Xtr=X[-ids,]; Xte=X[ids,]
    Utr=U[-ids,]; Ute=U[ids,]
    Ytr=Y[-ids,]; Yte=Y[ids,]
    betas=try(ridgeGLM2(Y=Ytr,X=Xtr,U=Utr, lambda = lambda, lambdaG = 0,
                        model = "surv", minSuccDiff = minSuccDiff,maxIter = maxIter))
    if(is(betas,"try-error")){
      lpnew=rep(0,nrow(Yte))} else {
        lpnew = cbind(Ute,Xte) %*% betas
      }
    RespAll <- rbind(RespAll,Yte)
    LPAll <- c(LPAll,as.numeric(lpnew))
  }
  
  CVLoss = .loglikSurv(Y=Surv(time=RespAll[,1], event = RespAll[,2]),lp=LPAll)
  return(-CVLoss)
}





.kcvLMloss_PgeqN_UP_GP <- function(lambdas, Y, X, U, Ds, folds, loss){
  
  ## ---------------------------------------------------------------------
  ## Internal function yields the cross-validated loglikelihood of the 
  ## ridge regression estimator with a two-dimensional covariate layout. 
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## lambdas : penalty parameter vector
  ## Y       : response vector
  ## X       : design matrix multiplied by the eigenvector matrix of the
  ##           nonnegative definite matrix that specifies the structure 
  ##           of spatial fused ridge penalty.
  ## U       : design matrix of the unpenalized covariates
  ## 
  ## Ds      : nonnegative eigenvalues of the nonnegative definite matrix 
  ##           that specifies the structure of spatial fused ridge penalty.
  ## folds   : list with folds
  ## loss    : character, either 'loglik' of 'sos', specifying the loss
  ##           criterion to be used in the cross-validation
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # evaluate loss per fold
  cvLoss <- 0
  
  # calculation prior to cv-loop
  X       <- sweep(X, 2, sqrt(lambdas[1] + lambdas[2] * Ds), FUN="/")
  XXT     <- tcrossprod(X)
  
  for (k in 1:length(folds)){
    # evaluate the estimator of the unpenalized low-dimensional 
    # regression parameter on the training samples
    if (all(dim(U) > 0)){
      gHat <- solve(crossprod(U[-folds[[k]],,drop=FALSE], 
                              solve(XXT[-folds[[k]], -folds[[k]]] + 
                                      diag(nrow(XXT)-length(folds[[k]])), 
                                    U[-folds[[k]],,drop=FALSE])),
                    crossprod(U[-folds[[k]],,drop=FALSE], 
                              solve(XXT[-folds[[k]], -folds[[k]]] + 
                                      diag(nrow(XXT)-length(folds[[k]])), 
                                    Y[-folds[[k]]])))
      gHat <- as.numeric(gHat)
    }
    if (all(dim(U) > 0)){
      Yhat   <- as.numeric(crossprod(t(XXT[folds[[k]], -folds[[k]]]),
                                     solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
                                             XXT[-folds[[k]], -folds[[k]]], 
                                           Y[-folds[[k]]] - as.numeric(crossprod(t(U[-folds[[k]],,drop=FALSE]), gHat))))) 
      + as.numeric(crossprod(t(U[-folds[[k]],,drop=FALSE]),gHat))
      
    } else {
      Yhat   <- 
        tcrossprod(solve(diag(rep(1, nrow(XXT)-length(folds[[k]]))) + 
                           XXT[-folds[[k]], -folds[[k]]], 
                         Y[-folds[[k]]]), 
                   XXT[folds[[k]], -folds[[k]], drop=FALSE])
    }
    
    
    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
        (log(2 * pi * sum((Y[folds[[k]]] - Yhat)^2) / 
               length(folds[[k]])) + 1) / 2
    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - Yhat)^2)
    }
  }
  
  # average over the folds
  return(cvLoss / length(folds))
}

.kcvLMloss_PleqN_UP_GPandNoGP <- function(lambdas, Y, X, U, Dg, folds, loss){
  
  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated log-likelihood of the ridge
  ## regression estimator  of the linear regression
  ## model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # evaluate loss per fold
  cvLoss <- 0
  for (k in 1:length(folds)){
    # evaluate regression estimate
    bHat <- ridgeGLM(Y[-folds[[k]]], X[-folds[[k]],,drop=FALSE],
                     U[-folds[[k]],,drop=FALSE], Dg=Dg,
                     lambda =lambdas[1], lambdaG=lambdas[2],
                     model="linear")
    
    # evaluate linear predictor		                 	
    Yhat <- X[folds[[k]],,drop=FALSE] %*% bHat[-c(1:ncol(U))] +
      U[folds[[k]],,drop=FALSE] %*% bHat[c(1:ncol(U))]
    
    if (loss == "loglik"){
      # evaluate the loglikelihood on the left-out samples
      cvLoss <- cvLoss + length(Y[folds[[k]]]) * 
        (log(2 * pi * sum((Y[folds[[k]]] - as.numeric(Yhat))^2) / 
               length(folds[[k]])) + 1) / 2
    }
    if (loss == "sos"){
      # evaluate the sum-of-sqaures on the left-out samples
      cvLoss <- cvLoss + sum((Y[folds[[k]]] - as.numeric(Yhat))^2)
    }
  }
  
  # average over the folds
  return(cvLoss / length(folds))
}

.kcvBLMloss_UP_GP <- function(lambdas, Y, X, U, Dg, folds, 
                              minSuccDiff, maxIter){
  
  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  ## 
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------
  
  # initiate
  cvLoss  <- 0
  lambda  <- exp(lambdas[1]) # positivity constraint
  lambdaG <- exp(lambdas[2]) # positivity constraint
  
  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){
    NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
    Dg <- tcrossprod(Matrix::kronecker(Diagonal(pX,1),
                                       solve(lambdaG*.PenMatr(NumNodes,1)+Diagonal(NumNodes,lambda)))
                     ,X)
    #Dg <- as.matrix(Dg)
    
    XDXT        <- crossprod(t(X), Dg)
    diagXDXTorg <- diag(XDXT)
  }
  if (ncol(X) < nrow(X)){ 
    Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
  }
  
  
  for (k in 1:length(folds)){
    
    
    
    
    # initialization
    penalty    <-  0
    lpOld <- rep(log(mean(Y[-folds[[k]]])/(1-mean(Y[-folds[[k]]]))), length(Y[-folds[[k]]]))
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    loglikPrev <- .loglikBLMlp(Y = Y[-folds[[k]]], lp = lpOld) #initial likelihood, penalty equals zero
    
    
    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lpOld))
      W0    <- Ypred * (1 - Ypred)
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <- 
          .Machine$double.eps
      }
      if (ncol(X) >= nrow(X)){ 
        # adjusted response
        Z <- lpOld + (Y[-folds[[k]]] - Ypred)/W0
        
        # evaluate unpenalized low-dim regression estimator
        diag(XDXT)[-folds[[k]]] <- diag(XDXT)[-folds[[k]]] + 1/W0
        slh        <- solve(XDXT[-folds[[k]], -folds[[k]]],Z) 
        
        gHat       <- solve(crossprod(U[-folds[[k]],], 
                                      solve(XDXT[-folds[[k]], -folds[[k]]], 
                                            U[-folds[[k]],])),
                            crossprod(U[-folds[[k]],], slh))
        Ug         <- tcrossprod(U, t(gHat))
        
        # evaluate linear predictor without evaluating the estimator
        slh        <- solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE], 
                            Z - Ug[-folds[[k]],])
        diag(XDXT)[-folds[[k]]] <- diagXDXTorg[-folds[[k]]]			
        lpAll      <- crossprod(XDXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- 0.5 * sum(lpAll[-folds[[k]]] * slh)
        lpAll      <- lpAll + Ug
      }
      
      if (ncol(X) < nrow(X)){ 		
        # adjusted response
        Z <- W0*lpOld + (Y[-folds[[k]]] - Ypred) 
        
        
        # evaluate subexpressions of the estimator
        XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 
                                 1, sqrt(W0), FUN="*")) 
        tUTX  <- crossprod(X[-folds[[k]], , drop=FALSE], 
                           sweep(U[-folds[[k]], , drop=FALSE], 
                                 1, W0, FUN="*"))
        XTZ   <- crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <- crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <- crossprod(sweep(U[-folds[[k]], , drop=FALSE], 
                                 1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
                                            crossprod(X[-folds[[k]], , drop=FALSE], 
                                                      sweep(U[-folds[[k]], , drop=FALSE], 
                                                            1, W0, FUN="*")))))))
        Ug   <- as.numeric(tcrossprod(gHat, U))
        
        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld + 
                      (Y[-folds[[k]]] - Ypred) - W0*Ug[-folds[[k]]]))
        bHat    <- qr.solve(XTXpD, XWZpT, 
                            tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
        penalty <- 0.5 * sum(crossprod(Dg, bHat) * 
                               (bHat))
      }
      
      # split linear predictor by fold
      lpOld      <- lpAll[-folds[[k]]]
      lpNew      <- lpAll[ folds[[k]]]
      
      # compute new likelihood
      loglik <- .loglikBLMlp(Y[-folds[[k]]], lpOld) - penalty
      
      # step-halving for stability
      if (loglik < loglikPrev){
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }
      
      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        #print("Singluar Fit")
        lpNew=rep(0,nrow(Y[folds[[k]]]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        #print("Succesful Fit")
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lpAll
      }
    }
    cvLoss <- cvLoss + .loglikBLMlp(Y[folds[[k]]], lpNew)
  }
  # average over the folds
  return(-cvLoss / length(folds))
}

.kcvSurvloss_UP_GP <- function(lambdas, Y, X, U, Dg, folds, 
                               minSuccDiff, maxIter){
  
  ## ---------------------------------------------------------------------
  ## Function yields the cross-validated loglikelihood of the ridge
  ## regression estimator of the logistic regression
  ## model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## folds   : list with folds
  ## 
  ## ...     : other arguments passed on to the 'ridgeGLM'-function
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the cross-validated loss of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # initiate
  RespTot <- c()
  LPTot <- c()
  
  lambda  <- exp(lambdas[1]) #positivity constraint
  lambdaG <- exp(lambdas[2]) #positivity constraint
  
  # evaluate subexpressions of the estimator
  if (ncol(X) >= nrow(X)){
    NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
    Dg <- tcrossprod(Matrix::kronecker(Diagonal(pX,1),
                                       solve(lambdaG*.PenMatr(NumNodes,1)+Diagonal(NumNodes,lambda)))
                     ,X)
    #Dg <- as.matrix(Dg)
    
    XDXT        <- crossprod(t(X), Dg)
    diagXDXTorg <- diag(XDXT)
  }
  if (ncol(X) < nrow(X)){ 
    Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
  }
  
  
  for (k in 1:length(folds)){
    
    # initialization
    penalty    <-  0
    lpPrev <- rep(0, length(Y[,2]))
    lpOld <- rep(0, length(Y[-folds[[k]],2]))
    loglikPrev <- .loglikSurv(Y=Y[-folds[[k]],], lp=lpOld)-penalty
    Yev = Y[-folds[[k]],2]
    
    
    for (iter in 1:maxIter){
      
      # calculate the weights
      H0 <- .breslow(Surv(time = Y[-folds[[k]],1], event = Y[-folds[[k]],2]),lpOld)[,2]
      Ypred <- as.numeric(H0 * exp(lpOld)) + 10^(-15)
      W0 <- Ypred 
      if (min(W0) <= 2*.Machine$double.eps){
        W0[which(W0 < 2*.Machine$double.eps)] <- 
          2*.Machine$double.eps
      }
      
      if (ncol(X) >= nrow(X)){ 
        # adjusted response
        Z <- lpOld + (Yev - Ypred)/W0
        
        # evaluate unpenalized low-dim regression estimator
        diag(XDXT)[-folds[[k]]] <- diag(XDXT)[-folds[[k]]] + 1/W0
        slh        <- solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE],Z) 
        
        gHat       <- solve(crossprod(U[-folds[[k]], , drop=FALSE], 
                                      solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE], 
                                            U[-folds[[k]], , drop=FALSE])),
                            crossprod(U[-folds[[k]], , drop=FALSE], slh))
        Ug         <- tcrossprod(U, t(gHat))
        
        # evaluate linear predictor without evaluating the estimator
        slh        <- solve(XDXT[-folds[[k]], -folds[[k]], drop=FALSE], 
                            Z - Ug[-folds[[k]],])
        diag(XDXT)[-folds[[k]]] <- diagXDXTorg[-folds[[k]]]			
        lpAll      <- crossprod(XDXT[-folds[[k]], ,drop=FALSE], slh)
        penalty    <- 0.5 * sum(lpAll[-folds[[k]]] * slh)
        lpAll      <- lpAll + Ug
      }
      
      if (ncol(X) < nrow(X)){ 		
        # adjusted response
        Z <- W0*lpOld + (Yev - Ypred) 
        
        
        # evaluate subexpressions of the estimator
        XTXpD <- crossprod(sweep(X[-folds[[k]], , drop=FALSE], 
                                 1, sqrt(W0), FUN="*")) 
        tUTX  <- crossprod(X[-folds[[k]], , drop=FALSE], 
                           sweep(U[-folds[[k]], , drop=FALSE], 
                                 1, W0, FUN="*"))
        XTZ   <- crossprod(X[-folds[[k]], , drop=FALSE], Z)
        UTZ   <- crossprod(U[-folds[[k]], , drop=FALSE], Z)
        UTU   <- crossprod(sweep(U[-folds[[k]], , drop=FALSE], 
                                 1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
                                            crossprod(X[-folds[[k]], , drop=FALSE], 
                                                      sweep(U[-folds[[k]], , drop=FALSE], 
                                                            1, W0, FUN="*")))))))
        Ug   <- as.numeric(tcrossprod(gHat, U))
        
        # evaluate linear predictor
        XWZpT <- as.numeric(
          crossprod(X[-folds[[k]], , drop=FALSE], W0*lpOld + 
                      (Yev - Ypred) - W0*Ug[-folds[[k]]]))
        bHat    <- qr.solve(XTXpD, XWZpT, 
                            tol=.Machine$double.eps)
        lpAll   <- as.numeric(tcrossprod(X, t(bHat))) + Ug
        penalty <- 0.5 * sum(crossprod(Dg, bHat) * 
                               (bHat))
      }
      
      # split linear predictor by fold
      lpOld      <- lpAll[-folds[[k]]]
      lpNew      <- lpAll[ folds[[k]]]
      
      #compute likelihood
      loglik <- .loglikSurv(Surv(time=Y[-folds[[k]],1],event = Y[-folds[[k]],2]), lpOld) - penalty
      
      # step-halving for stability
      if (loglik < loglikPrev){ 
        lpOld <- 0.5*lpOld + 0.5*lpPrev[-folds[[k]]]
        lpNew <- 0.5*lpNew + 0.5*lpPrev[folds[[k]]]
      }
      
      # assess convergence
      if (is.infinite(loglik) | is.na(loglik) | is.nan(loglik)){
        print("Singluar Fit")
        lpNew=rep(0,nrow(Y[folds[[k]],]))
        break
      } else if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        #print("Succesful Fit")
        break
      } else {
        loglikPrev <- loglik
        lpPrev <- lpAll
      }
      if(iter==maxIter){
        print("Not converged yet, please increase maxIter")
        lpNew=rep(0,nrow(Y[folds[[k]],]))
      }
    }
    
    RespTot <- base::append(RespTot,Y[folds[[k]],])
    LPTot <- c(LPTot,as.numeric(lpNew))
    
  }
  CVLoss = .loglikSurv(Y=Surv(time=RespTot[,1],event = RespTot[,2]),lp=LPTot)
  # average over the folds
  return(-CVLoss)
  
}
.kcvSurvloss_UP_GP1<- function(lambdas, Y, X, U, Dg, folds, 
                               minSuccDiff, maxIter){
  
  RespAll <- c()
  LPAll <- c()
  lambda  <- lambdas[1]
  lambdaG <- lambdas[2]
  
  for (k in 1:length(folds)){
    ids=folds[[k]]
    Xtr=X[-ids,]; Xte=X[ids,]
    Utr=U[-ids,]; Ute=U[ids,]
    Ytr=Y[-ids,]; Yte=Y[ids,]
    betas=try(ridgeGLM2(Y=Ytr,X=Xtr,U=Utr, lambda = lambda, lambdaG = lambdaG,
                        Dg = Dg, model = "surv", minSuccDiff = minSuccDiff,maxIter = maxIter))
    if(is(betas,"try-error")){
      lpnew=rep(0,nrow(Yte))} else {
        lpnew = cbind(Ute,Xte) %*% betas
      }
    RespAll <- rbind(RespAll,Yte)
    LPAll <- c(LPAll,as.numeric(lpnew))
  }
  
  CVLoss = .loglikSurv(Y=Surv(time=RespAll[,1],event = RespAll[,2]),lp=LPAll)
  return(-CVLoss)
}





.ridgeLM <- function(Y, X, U, lambda, lambdaG, Dg){
  ## ---------------------------------------------------------------------
  ## Function that evaluates ridge regression estimator with a regular and
  ## generalized penalty. The fused penalty is specified by the matrix Dg.
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix multiplied by the eigenvector matrix of the
  ##           nonnegative definite matrix that specifies the structure 
  ##           of spatial fused ridge penalty.
  ## lambda  : numeric, the regular ridge penalty parameter
  ## lambdaG : numeric, the fused ridge penalty parameter
  ## Dg      : nonnegative definite matrix that specifies the structure 
  ##           of the generalized ridge penalty.
  
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the generalized ridge regression estimate.
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen, modified by JM Goedhart
  ## ---------------------------------------------------------------------
  
  
  if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) > 0)){
    # no generalized ridge penalization + unpenalized covariates
    
    if (nrow(X) >= ncol(X)){
      # efficient evaluation for n >= p
      
      # evaluate subexpressions of the estimator
      tUTX        <- crossprod(X, U)
      XTXpI       <- crossprod(X) 
      XTY         <- crossprod(X, Y) 
      UTY         <- crossprod(U, Y) 
      UTU         <- crossprod(U)
      diag(XTXpI) <- diag(XTXpI) + lambda
      
      # evaluate unpenalized low-dim regression estimator
      gHat <- solve(UTU - crossprod(solve(XTXpI, tUTX), tUTX), 
                    (UTY - as.numeric(crossprod(solve(XTXpI, XTY), tUTX))))
      
      # evaluate penalized high-dim regression estimator
      bHat <- solve(XTXpI, XTY - crossprod(t(tUTX), gHat))
    }
    
    if (nrow(X) < ncol(X)){
      # efficient evaluation for n < p
      
      # evaluate subexpressions of the estimator
      XXT       <- tcrossprod(X) / lambda
      diag(XXT) <- diag(XXT) + 1
      Y         <- Y 
      
      # evaluate unpenalized low-dim regression estimator
      gHat <- solve(crossprod(U, solve(XXT, U)),
                    crossprod(U, solve(XXT, Y)))
      
      # evaluate penalized high-dim regression estimator
      Y    <- Y - crossprod(t(U), gHat) 
      bHat <- crossprod(X, solve(XXT, Y))/lambda
    }
  }
  
  if (((max(abs(Dg)) > 0) && lambdaG > 0) && (ncol(U) > 0)){
    if (nrow(X) >= ncol(X)){
      # efficient evaluation for n >= p
      
      # evaluate subexpressions of the estimator
      Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
      XTX  <- crossprod(X)
      tUTX <- crossprod(X, U)
      XTY  <- crossprod(X, Y)
      UTY  <- crossprod(U, Y)
      UTU  <- crossprod(U)
      
      # evaluate unpenalized low-dim regression estimator
      gHat <- solve(UTU -  crossprod(solve(Dg + XTX, tUTX), tUTX), 
                    (UTY - crossprod(solve(Dg + XTX, XTY), tUTX)[1,]))
      
      # evaluate penalized high-dim regression estimator
      bHat <- solve(XTX + Dg, XTY - crossprod(t(tUTX), gHat))
    }
    
    if (nrow(X) < ncol(X)){
      # efficient evaluation for n < p
      
      # evaluate subexpressions of the estimator
      NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
      Dg <-  tcrossprod(Matrix::kronecker(Diagonal(pX,1),
                                          solve(lambdaG*.PenMatr(NumNodes,1)+Diagonal(NumNodes,lambda)))
                        ,X)
      #Dg <- as.matrix(Dg)
      
      
      XDXTpI <- crossprod(t(X), Dg) + Diagonal(nrow(X),1)
      
      # evaluate unpenalized low-dim regression estimator
      gHat <- solve(crossprod(U, solve(XDXTpI, U)),
                    crossprod(U, solve(XDXTpI, Y)))
      
      # evaluate penalized high-dim regression estimator
      Y    <- Y - crossprod(t(U), gHat) 
      bHat <- tcrossprod(Dg, t(solve(XDXTpI, Y)))[,1]
    }
  }
  return(c(as.numeric(gHat), as.numeric(bHat)))
}


.ridgeBLM <- function(Y, X, U, lambda, lambdaG, Dg, minSuccDiff, maxIter){
  
  if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) > 0)){
    # no generalized ridge penalization + unpenalized covariates
    
    # initiate
    
    if (ncol(X) >= nrow(X)){ 
      XXT <- tcrossprod(X) / lambda
    }
    if (ncol(X) >= nrow(X)){ 
      tUTX <- crossprod(X, U)
    }
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    lp      <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    loglikPrev <- .loglikBLMlp(Y, lp)
    
    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lp))
      W0    <- as.numeric(Ypred * (1 - Ypred))
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <- 
          .Machine$double.eps
      }
      
      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){ 
        # obtain part one of the estimator
        Z <- lp + (Y - Ypred)/W0
        
        # now obtain the IRWLS update efficiently
        diag(XXT) <- diag(XXT) + 1/W0
        slh       <- solve(XXT, Z)
        gHat      <- solve(crossprod(U, solve(XXT, U)),
                           crossprod(U, slh))
        slh       <- solve(XXT, Z - U %*% gHat)				
        diag(XXT) <- diag(XXT) - 1/W0
        lp        <- crossprod(XXT, slh)
        penalty   <- sum(lp * slh) / 2
        lp        <- lp + U %*% gHat
        loglik    <- .loglikS(Y, lp) - penalty
      }
      if (ncol(X) < nrow(X)){ 
        # adjusted response
        Z <- W0*lp + (Y - Ypred)
        
        # evaluate subexpressions of the estimator
        XTXpI       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
        tUTX        <- crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ         <- crossprod(X, Z)
        UTZ         <- crossprod(U, Z)
        UTU         <- crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        diag(XTXpI) <- diag(XTXpI) + lambda
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpI, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpI, XTZ), 
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))
        
        # obtain part one of the estimator
        XWZpT <- as.numeric(
          crossprod(X, W0*lp + (Y - Ypred) - 
                      W0*as.numeric(tcrossprod(gHat, U))))
        
        # now obtain the IRWLS update efficiently
        bHat <- solve(XTXpI, XWZpT)
        lp   <- as.numeric(tcrossprod(X, t(bHat)) + 
                             as.numeric(tcrossprod(gHat, U)))
        
        # evaluate the loglikelihood
        loglik <- .loglikBLMlp(Y, lp) -
          0.5 * lambda * sum((bHat)^2)
      }
      print(loglik)
      #step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev 
      }
      
      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)){stop("convergence error, please increase penalty")}
      
      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        print(paste("IRLS converged at iteration ",iter))
        break 
      } else {
        loglikPrev <- loglik
        lpPrev <-lp
      }
      if(iter==maxIter){print("Not converged yet, please increase maxIter")}
    }
    if (ncol(X) >= nrow(X)){ 
      bHat <- crossprod(X, slh) / lambda
    }
  }
  
  if (((max(abs(Dg)) != 0) && lambdaG != 0) && (ncol(U) > 0)){
    
    if (ncol(X) >= nrow(X)){ 
      # evaluate subexpressions of the estimator
      NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
      Dg <- tcrossprod(Matrix::kronecker(Diagonal(pX,1),
                                         solve(lambdaG*.PenMatr(NumNodes,1)+Diagonal(NumNodes,lambda)))
                       ,X)
      #Dg <- as.matrix(Dg)
      
      XDXT        <- crossprod(t(X), Dg)
      diagXDXTorg <- diag(XDXT)
    } else { 
      # evaluate subexpressions of the estimator
      Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
    }
    
    # initiate
    lpPrev <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    lp      <- rep(log(mean(Y)/(1-mean(Y))), length(Y))
    loglikPrev <- .loglikBLMlp(Y, lp)
    
    
    for (iter in 1:maxIter){
      # calculate the weights
      Ypred <- 1 / (1 + exp(-lp))
      W0    <- as.numeric(Ypred * (1 - Ypred))
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <- 
          .Machine$double.eps
      }
      
      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){ 
        # obtain part one of the estimator
        Z <- lp + (Y - Ypred)/W0
        
        # now obtain the IRWLS update efficiently
        diag(XDXT) <- diag(XDXT) + 1/W0
        slh        <- qr.solve(XDXT, Z, 
                               tol=.Machine$double.eps)
        gHat       <- solve(crossprod(U, solve(XDXT, U)),
                            crossprod(U, slh))
        slh        <- qr.solve(XDXT, Z - U %*% gHat,
                               tol=.Machine$double.eps)
        diag(XDXT) <- diagXDXTorg
        lp         <- crossprod(XDXT, slh)
        penalty    <- sum(lp * slh) / 2
        lp         <- lp + U %*% gHat
        loglik    <- .loglikBLMlp(Y, lp) - penalty
      }
      if (ncol(X) < nrow(X)){ 
        # adjusted response
        Z <- W0*lp + (Y - Ypred)
        
        # evaluate subexpressions of the estimator
        XTXpD <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
        tUTX  <- crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ   <- crossprod(X, Z)
        UTZ   <- crossprod(U, Z)
        UTU   <- crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))
        
        # obtain part one of the estimator
        XWZpT <-  
          crossprod(X, W0*lp + (Y - Ypred) - 
                      W0*tcrossprod(gHat, U)[1,])[,1]
        
        # now obtain the IRWLS update efficiently
        bHat    <- qr.solve(XTXpD, XWZpT, 
                            tol=.Machine$double.eps)
        penalty <- sum(crossprod(Dg, bHat) * 
                         (bHat))/2
        lp      <- as.numeric(tcrossprod(X, t(bHat)) + 
                                as.numeric(tcrossprod(gHat, U)))
        
        # evaluate the loglikelihood
        loglik    <- .loglikBLMlp(Y, lp) - penalty
      }
      
      # step-halving for stability
      print(loglik)
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev #step-halving for first iteration to guide search
      }
      
      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)){stop("convergence error, please increase penalty")}
      
      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        print(paste("IRLS converged at iteration ",iter))
        break 
      } else {
        loglikPrev <- loglik
        lpPrev <- lp
      }
      if(iter==maxIter){print("Not converged yet, please increase maxIter")}
    }
    if (ncol(X) >= nrow(X)){ 
      bHat <- as.numeric(tcrossprod(as.numeric(slh), Dg))
    }
  }
  return(c(as.numeric(gHat), as.numeric(bHat)))
}

.ridgeSurv <- function(Y,X,U,lambda,lambdaG, Dg, minSuccDiff,maxIter){
  
  
  if (((max(abs(Dg)) == 0) | lambdaG == 0) && (ncol(U) > 0)){
    # no generalized ridge penalization + unpenalized covariates
    
    # initiate
    
    if (ncol(X) >= nrow(X)){ 
      XXT <- tcrossprod(X) /lambda
    }
    if (ncol(X) >= nrow(X)){ 
      tUTX <- crossprod(X, U)
    }
    
    #initiate
    lpPrev <- rep(0, length(Y))
    lp <- rep(0, length(Y))
    loglikPrev <- .loglikSurv(Y, lp) #penalty is zero for initial betas
    print(loglikPrev)
    Yev = Y[,2] # event part of response
    
    for (iter in 1:maxIter){
      # calculate the weights
      
      H0 <- .breslow(Y,lp)[,2]
      Ypred <- as.numeric(H0 * exp(lp))
      W0 <- Ypred 
      
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <- 
          2*.Machine$double.eps}
      
      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){ 
        # obtain part one of the estimator
        Z <- lp + (Yev - Ypred)/W0
        
        # now obtain the IRWLS update efficiently
        diag(XXT) <- diag(XXT) + 1/W0
        slh       <- solve(XXT, Z)
        gHat      <- solve(crossprod(U, solve(XXT, U)),
                           crossprod(U, slh))
        
        
        slh       <- solve(XXT, Z - U %*% gHat)	
        #Z1 <- W0*lp + Yev-Ypred - diag(W0)%*% U %*% gHat
        diag(XXT) <- diag(XXT) - 1/W0
        #Hatmat <- XXT-XXT%*%solve(XXT+diag(1/W0),XXT)
        #lp <- Hatmat %*% Z1
        
        lp        <- crossprod(XXT, slh)
        penalty   <- sum(lp * slh) / 2
        lp        <- lp + U %*% gHat
      }
      if (ncol(X) < nrow(X)){ 
        # adjusted response
        Z <- W0*lp + (Yev - Ypred)
        
        # evaluate subexpressions of the estimator
        XTXpI       <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
        tUTX        <- crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ         <- crossprod(X, Z)
        UTZ         <- crossprod(U, Z)
        UTU         <- crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        diag(XTXpI) <- diag(XTXpI) + lambda
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpI, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpI, XTZ), 
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))
        
        # obtain part one of the estimator
        XWZpT <- as.numeric(
          crossprod(X, W0*lp + (Yev - Ypred) - 
                      W0*as.numeric(tcrossprod(gHat, U))))
        
        # now obtain the IRWLS update efficiently
        bHat <- qr.solve(XTXpI, XWZpT, tol = .Machine$double.eps)
        lp   <- as.numeric(tcrossprod(X, t(bHat)) + 
                             as.numeric(tcrossprod(gHat, U)))
        penalty <- 0.5 * lambda * sum((bHat)^2)
      }
      
      loglik    <- .loglikSurv(Y, lp) - penalty
      print(loglik)
      # step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev 
      }
      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)){stop("convergence error, please increase penalty")}
      
      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        print(paste("IRLS converged at iteration ",iter))
        break 
      } else {
        loglikPrev <- loglik
        lpPrev <- lp
      }
      if(iter==maxIter){
        print("Not converged yet, please increase maxIter")}
      
    }
    
    if (ncol(X) >= nrow(X)){ 
      bHat <- crossprod(X, slh) / lambda
    }
  }
  
  if (((max(abs(Dg)) != 0) && lambdaG != 0) && (ncol(U) > 0)){
    
    
    if (ncol(X) >= nrow(X)){ 
      # evaluate subexpressions of the estimator
      NumNodes = length(which(Dg[1,]!=0)); pX=ncol(Dg)/NumNodes
      Dg <- tcrossprod(Matrix::kronecker(Diagonal(pX,1),
                                         solve(lambdaG*.PenMatr(NumNodes,1)+Diagonal(NumNodes,lambda)))
                       ,X)
      #Dg <- as.matrix(Dg)
      
      XDXT        <- crossprod(t(X), Dg)
      diagXDXTorg <- diag(XDXT)
    } else { 
      # evaluate subexpressions of the estimator
      Dg   <- diag(rep(lambda, ncol(X))) + lambdaG * Dg
    }
    
    #initiate
    lpPrev <- rep(0, length(Y))
    lp <- rep(0, length(Y))
    loglikPrev <- .loglikSurv(Y, lp) #penalty is zero for initial betas
    Yev=Y[,2]
    
    for (iter in 1:maxIter){
      # calculate the weights
      
      H0 <- .breslow(Y,lp)[,2]
      Ypred <- as.numeric(H0 * exp(lp))
      W0 <- Ypred 
      
      if (min(W0) <= .Machine$double.eps){
        W0[which(W0 < .Machine$double.eps)] <- 
          2*.Machine$double.eps}
      
      # distinguish between the high- and low-dimensional cases
      if (ncol(X) >= nrow(X)){ 
        # obtain part one of the estimator
        Z <- lp + (Yev - Ypred)/W0
        
        # now obtain the IRWLS update efficiently
        diag(XDXT) <- diag(XDXT) + 1/W0
        slh        <- qr.solve(XDXT, Z, 
                               tol=.Machine$double.eps)
        gHat       <- solve(crossprod(U, qr.solve(XDXT, U,tol = .Machine$double.eps)),
                            crossprod(U, slh))
        slh        <- qr.solve(XDXT, Z - U %*% gHat,
                               tol=.Machine$double.eps)
        diag(XDXT) <- diagXDXTorg
        lp         <- crossprod(XDXT, slh)
        penalty    <- sum(lp * slh) / 2
        lp         <- lp + U %*% gHat
      }
      if (ncol(X) < nrow(X)){ 
        # adjusted response
        Z <- W0*lp + (Yev - Ypred)
        
        # evaluate subexpressions of the estimator
        XTXpD <- crossprod(sweep(X, 1, sqrt(W0), FUN="*")) 
        tUTX  <- crossprod(X, sweep(U, 1, W0, FUN="*"))
        XTZ   <- crossprod(X, Z)
        UTZ   <- crossprod(U, Z)
        UTU   <- crossprod(sweep(U, 1, sqrt(W0), FUN="*"))
        XTXpD <- XTXpD + Dg
        
        # evaluate unpenalized low-dim regression estimator
        gHat <- as.numeric(
          solve(UTU - crossprod(solve(XTXpD, tUTX), tUTX), 
                (UTZ - as.numeric(crossprod(solve(XTXpD, XTZ), 
                                            crossprod(X, sweep(U, 1, W0, FUN="*")))))))
        
        # obtain part one of the estimator
        XWZpT <-  
          crossprod(X, W0*lp + (Yev - Ypred) - 
                      W0*tcrossprod(gHat, U)[1,])[,1]
        
        # now obtain the IRWLS update efficiently
        bHat    <- qr.solve(XTXpD, XWZpT, 
                            tol=.Machine$double.eps)
        penalty <- sum(crossprod(Dg, bHat) * 
                         (bHat))/2
        lp      <- as.numeric(tcrossprod(X, t(bHat)) + 
                                as.numeric(tcrossprod(gHat, U)))
      }
      # compute new likelihood
      loglik    <- .loglikSurv(Y, lp) - penalty
      print(loglik)
      # step-halving for stability
      if (loglik < loglikPrev){
        lp <- 0.5*lp + 0.5*lpPrev 
      }
      
      # assess convergence
      if(is.nan(loglik) | is.infinite(loglik) | is.na(loglik)){stop("convergence error, please increase penalty")}
      
      if (abs((loglik - loglikPrev)/loglikPrev) < minSuccDiff){
        print(paste("IRLS converged at iteration ",iter))
        break 
      } else {
        loglikPrev <- loglik
        lpPrev <- lp
      }
      if(iter==maxIter){
        print("Not converged yet, please increase maxIter")}
    }
    
    if (ncol(X) >= nrow(X)){ 
      bHat <- as.numeric(tcrossprod(as.numeric(slh), Dg))
    }
  }
  return(c(as.numeric(gHat), as.numeric(bHat)))
}



###### Auxiliary functions ######
#################################

.loglikBLM <- function(Y, X, betas){
  
  ## ---------------------------------------------------------------------
  ## Function calculates the loglikelihood of the logistic regression model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## betas   : regression parameter
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the loglikelihood of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # evaluate the linear predictor
  lp <- as.numeric(tcrossprod(betas, X))
  
  # evaluate the loglikelihood
  loglik1 <- exp(Y * lp) / (1 + exp(lp))
  loglik2 <- exp(((Y-1) * lp)) / (1 + exp(-lp))
  loglik1[!is.finite(log(loglik1))] <- NA
  loglik2[!is.finite(log(loglik2))] <- NA
  loglik  <- sum(log(apply(cbind(loglik1, loglik2), 1, mean, na.rm=TRUE)))
  return(loglik)
}

.loglikBLMlp <- function(Y, lp){
  
  ## ---------------------------------------------------------------------
  ## Function calculates the loglikelihood of the logistic regression model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## lp      : linear predictor
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the loglikelihood of the logistic regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  # evaluate the loglikelihood
  loglik1 <- exp(Y * lp) / (1 + exp(lp))
  loglik2 <- exp(((Y-1) * lp)) / (1 + exp(-lp))
  loglik1[!is.finite(log(loglik1))] <- NA
  loglik2[!is.finite(log(loglik2))] <- NA
  loglik  <- sum(log(apply(cbind(loglik1, loglik2), 1, mean, na.rm=TRUE)))
  return(loglik)
}

.loglikLM <- function(Y, X, betas){
  
  ## ---------------------------------------------------------------------
  ## Function calculates the loglikelihood of the linear regression model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## betas   : regression parameter
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the loglikelihood of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  return(-length(Y) * (log(2 * pi * sum((Y - as.numeric(tcrossprod(as.numeric(betas), X)))^2) / length(Y)) + 1) / 2)
}

.breslow <- function(Y,lp){ #checked to coincide with basehaz from penalized (which coincides with survfit)
  #Y survival response; lp: linear predictor (length n)
  #Returns hazard and cumulative hazard at all timepoints
  # Author: Mark A van de Wiel
  if(class(Y)!="Surv"){stop("Survival Model, please specify Y as survival object")}
  ord <- order(Y)
  invord <- order(ord) #to put back in original order
  Ys <- Y[ord]
  di <- Ys[,2] #event
  ti <- Ys[,1] #time
  lp <- lp[ord]
  htsort <- di/rev(cumsum(exp(rev(lp))))
  if (min(htsort) <= .Machine$double.eps){
    htsort[which(htsort < .Machine$double.eps)] <- 
      2*.Machine$double.eps}
  #ht <- htsort
  #Ht <- cumsum(htsort)
  ht <- htsort[invord]
  Ht <- cumsum(htsort)[invord]
  return(data.frame(ht=ht,Ht=Ht,time=ti[invord]))
}

.loglikSurv <- function(Y,lp){
  #Author: Mark A. van de Wiel
  if(class(Y)!="Surv"){stop("Survival Model, please specify Y as survival object")}
  hazards <- .breslow(Y,lp);
  di <- Y[,2]
  ht <- hazards[,1]
  Ht <- hazards[,2]
  #print(ht)
  #print(Ht)
  thescore <- sum(-Ht*exp(lp))
  di1 <- which(di==1)
  if(length(di1)>0) thescore <- thescore + sum(di[di1]*(log(ht[di1])+lp[di1]))
  return(thescore)
}

.CIndexSurv <- function(Y,lp){
  # Author: Mark A van de Wiel
  if(class(Y)!="Surv"){stop("Survival Model, please specify Y as survival object")}
  lpmin <- -lp
  thescore <- concordance(response ~ lpmin)$concordance
  if(is.na(thescore)) thescore <- 0.5
  return(thescore)
}

LogLik <- function(Y,X,betas,model){
  if (model=="linear"){
    return(.loglikLM(Y,X,betas))
  }
  if (model=="logistic"){
    return(.loglikBLM(Y,X,betas))
  }
}

.sosLM <- function(Y, X, betas){
  
  ## ---------------------------------------------------------------------
  ## Function calculates the sum-of-squares of the linear regression model
  ## --------------------------------------------------------------------- 
  ## Arguments
  ## Y       : response vector
  ## X       : design matrix
  ## betas   : regression parameter
  ## ---------------------------------------------------------------------
  ## Value:
  ## A numeric, the sum-of-squares of the linear regression model
  ## ---------------------------------------------------------------------
  ## Authors : Wessel N. van Wieringen
  ## ---------------------------------------------------------------------
  
  return(sum((Y - as.numeric(tcrossprod(as.numeric(betas), X)))^2))
}



