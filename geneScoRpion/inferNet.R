inferNet <- function(Y, P, method="lsco", zetas=linspace(0,1,100), zetas.auto=FALSE){

  handleZeta <- function(A, zetas=linspace(0,1,100), zetas.auto2=zetas.auto){
    if(zetas.auto2)
    {
      zetas <- linspace(0,max(abs(A)),100)
    }
    A <- as.matrix(A)
    zetasMin <- min(abs(A[A>0]))-.Machine$double.eps#eps(1)
    zetasMax <- max(abs(A[A>0]))+.Machine$double.eps*10#10*eps(1)
    delta <- zetasMax - zetasMin
    zetas <- zetas*delta+zetasMin
    
    Alst <- list()
    for(i in 1:length(zetas)){
      Atmp <- A
      Atmp[abs(Atmp) <= zetas[i]] <- 0
      Alst[[i]] <- Atmp
    }
    return(Alst)
  }
  
  handleZeta2 <- function(A, zetas=linspace(0,1,100)){
    A <- as.matrix(A)
    zetasMin <- min(abs(A[A>0]))-.Machine$double.eps#eps(1)
    zetasMax <- max(abs(A[A>0]))+.Machine$double.eps*10#10*eps(1)
    delta <- zetasMax - zetasMin
    zetas <- zetas*delta + zetasMin
    
    Alst <- list()
    for(i in 1:length(zetas)){
      Atmp <- A
      Atmp[abs(Atmp) <= zetas[i]] <- 0
      Alst[[i]] <- Atmp
    }
    return(Alst)
  }
  
  ##### none #####
  if(method == "none"){
    A <- Y
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  ##### lsco #####
  if(method == "lsco"){
    library(pracma)
    pv <- pinv(as.matrix(Y))
    A <- -P %*% pinv(Y)
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  ##### lasso #####
  if(method == "lasso"){
  library(glmnet)
  
  A <- list()
  for(i in 1:dim(P)[1]){
    fit <- glmnet(t(Y), as.matrix(-P[i,]), lambda = zetas, alpha = 1)
  
  betas <- as.matrix(fit$beta)
  betas2 <- betas[,ncol(betas):1]
  A[[i]] <- betas2
  }
  Az <- list()
  for(i in 1:length(zetas)){
    Az[[i]] <- do.call(rbind, lapply(A, function(x) x[,i]))#as.matrix(sapply(do.call(rbind, lapply(A, function(x) x[,i])), as.numeric))  
  }
  return(Az)
  }
  
  ##### ridge #####
  if(method == "ridge"){
    library(glmnet)
    
    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      fit <- glmnet(t(Y), as.matrix(-P[i,]), alpha = 0)
      A[i,] <- as.numeric(as.matrix(coef(fit, s = eps(1))))[-1]
    }
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  ##### lsfit #####
  if(method == "lsfit"){

    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      fit <- lsfit(t(Y), as.matrix(-P[i,]), intercept=F)
      A[i,] <- as.numeric(fit$coefficients)
    }
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  ##### lm #####
  if(method == "lm"){
    
    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      fit <- lm(t(Y) ~ as.matrix(-P[i,]))
      A[i,] <- fit$coefficients[2,]
    }
    Az <- handleZeta(t(A), zetas) #flip it
    return(Az)
  }
  
  
  #### QR ####
  if(method == "qr"){
    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
    A[i,] <- qr.solve(t(Y),-P[i,])
    }
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  #### LOESS ####
  if(method == "loess"){
    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      #data <- data.frame(t(Y), decision = -P[i,])
      data <- data.frame(Y[,seq(i,i+100*dim(P)[1],dim(P)[1])[1:sum(P[1,])]], decision=-as.matrix(P[i,1:dim(P)[1]]))
      fit <- loess(decision~., data, span=0.1, drop.square = T)
      A[i,] <- fit$fitted
    }
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  
  if(method == "msir"){
    library(msir)
    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      x <- as.matrix(t(Y))
      y <- -as.matrix(as.numeric(P[i,]))
      fit <- msir(x,y,G=1)
      A[i,] <- as.numeric(fit$basis) #fit$basis  fit$evalues
    }
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  if(method == "sir"){
    library(dr)
    A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      x <- as.matrix(t(Y))
      y <- -as.matrix(as.numeric(P[i,]))
      fit <- dr::dr.compute(x,y,rep(1,dim(x)[1]), method = "save")
      A[i,] <- as.numeric(fit$result[[1]]$B) #as.numeric(fit$evalues) #fit$basis  fit$evalues
    }
    Az <- handleZeta(A, zetas)
    return(Az)
  }

  ## RANDOM FOREST
  
  if(method == "rf"){
    #library(randomForest)
    library(ranger)
    #A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
    #for(i in 1:dim(P)[1]){
      #i=1
      #data <- data.frame(Y, decision=as.matrix(P[i,1:dim(P)[1]]))
      #data <- data.frame(Y[,seq(i,i+100*dim(P)[1],dim(P)[1])[1:sum(P[1,])]], decision=as.matrix(P[i,1:dim(P)[1]]))
      #rf <- randomForest(decision~., data, proximity=TRUE)
      #A[i,] <- unname(rf$predicted)
    #}
    A <- matrix(0,ncol = dim(P)[1], nrow = dim(P)[1])
    for(i in 1:dim(P)[1]){
      #i=1
      data <- data.frame(t(Y), decision=-as.matrix(P[i,]))
      #data <- data.frame(Y, decision=-as.matrix(P[i,1:dim(P)[1]]))
      #data <- data.frame(Y[,seq(i,i+100*dim(P)[1],dim(P)[1])[1:sum(P[1,])]], decision=as.matrix(P[i,1:dim(P)[1]]))
      
      #rf <- ranger(decision~., data, case.weights=-(as.matrix(P[i,]))/2+0.2, importance = "impurity",
      #             classification = F, num.trees = 10000, splitrule = "variance", holdout = F, replace = T)
      
      rf <- ranger(decision~., data, importance = "impurity", probability = T,
                   classification = T, num.trees = 1000, splitrule = "extratrees", replace = F)
      

      A[i,] <- unname(rf$variable.importance)
    }
    #A[A==0] <- 1
    A <- scale(A)
    A <- t(-A)
    Az <- handleZeta(A, zetas)
    return(Az)
  }

  #### FOR 1 REP ####
  
  ##### genie3 #####
  if(method == "GENIE3"){
    library(GENIE3)
    
    if(is.null(rownames(Y))){
      rownames(Y) <- paste0(rep("G",dim(Y)[1]),1:dim(Y)[1])
    }
    A <- GENIE3(Y)
    #diag(A) <- 1
    A2 <- t(A)
    Az <- handleZeta(A2, zetas)
    return(Az)
  }
  
  ##### neural networks #####
  if(method == "lnnet"){
    library(nnet)
    l <- list()
    for(i in 1:dim(Y)[1]){ # 100 is how many replicates in total (can be increased)
      nnet.fit <- nnet(Y[,seq(i,i+100*dim(P)[1],dim(P)[1])[1:sum(P[1,])]], as.matrix(P[i,1:dim(P)[1]]),
                       size= 1, linout = T, skip = T, Hess = T)
      #nnet.fit <- nnet(t(Y) ~ as.matrix(-P[i,]), size=1, linout = T)
      l[[i]] <- nnet.fit$fitted.values
    }
    A <- do.call(cbind.data.frame, l)
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  ##### neural networks 2 #####
  if(method == "lnnet2"){
    library(neuralnet)
    
    l <- list()
    for(i in 1:dim(Y)[1]){
      data <- data.frame(Y[,seq(i,i+100*dim(P)[1],dim(P)[1])[1:sum(P[1,])]], decision=-as.matrix(P[i,1:dim(P)[1]]))
      nnet.fit2 <- neuralnet(decision ~., data, linear.output = T, hidden = c(dim(P)[1]), algorithm = 'rprop+')
      l[[i]] <- unlist(nnet.fit2$net.result)
      print(i)
    }
      A <- do.call(cbind.data.frame, l)
      colnames(A) <- NULL
      Az <- handleZeta(A, zetas, zetas.auto2=F)
      return(Az)
  }
  
  ##### SVM #####
  
  if(method == "svm"){
    library(e1071)
    
    l <- list()
    for(i in 1:dim(Y)[1]){
      data <- data.frame(Y[,seq(i,i+100*dim(P)[1],dim(P)[1])[1:sum(P[1,])]], decision=as.matrix(P[i,1:dim(P)[1]]))
      svmfit <- svm(decision ~., data, type="nu-regression", kernel="sigmoid")
      l[[i]] <- unlist(svmfit$decision.values)
    }
    A <- do.call(cbind.data.frame, l)
    colnames(A) <- NULL
    Az <- handleZeta(A, zetas)
    return(Az)
  }
  
  
}