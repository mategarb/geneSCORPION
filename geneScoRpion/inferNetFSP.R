inferNetFSP <- function(Y, P, fs, fsmeth, topf, permn, pval){
  #Y <- Y %>% as.data.frame
 # P <- P %>% as.data.frame
A <- matrix(0,ncol = dim(P)[1],nrow = dim(P)[1])
#outB <- Boruta(dec~.,data=Y4,doTrace=2,pValue=1,mcAdj=F)
#outMARS <- earth(dec~.,data=Y4)


if(fs){
  for(i in 1:dim(P)[1]){
    Y2 <- data.frame(t(Y), dec=-P[i,])
    colnames(Y2)[dim(Y2)[2]] <- "dec" 
    dd <- mRMR.data(data = Y2)
    # preselect only N informative genes
    
    fsres <- mRMR.classic(data = dd, target_indices = dim(Y2)[2], feature_count = topf)
    ids <- solutions(fsres) %>% unlist %>% unname
    # infer with lasso
    fit <- glmnet(Y2[,ids], as.matrix(-P[i,]), alpha = 1)
    tmpA <- as.matrix(coef(fit, s = eps(1)))[-1,]
    A[match(names(tmpA),colnames(Y2)),i] <- unname(tmpA)
    
    ## perm test
    tmpAS <- list()
    for (p in 1:permn){
      fitS <- glmnet(Y2[,ids], sample(as.matrix(-P[i,])), alpha = 1)
      tmpAS[[p]] <- as.matrix(coef(fitS, s = eps(1)))[-1,]
    }
    tmpAS2 <- do.call(rbind.data.frame, tmpAS)
    p <- c()
    for(k in 1:length(tmpA)){
      if(tmpA[k] <= 0){
        p[k] <- (length(which(tmpAS2[,k] <= as.numeric(tmpA[k])))/permn)+1/permn
      }else
      {
        p[k] <- (length(which(tmpAS2[,k] >= as.numeric(tmpA[k])))/permn)+1/permn
      }
    }
    A[match(names(tmpA)[which(p >= pval)],colnames(Y2)),i] <- 0
    progress(i,dim(P)[1], progress.bar = T)
    if(i == dim(P)[1]) message("Done!")
  }
  diag(A) <- 0
  A <- A %>% as.data.frame
  colnames(A) <- colnames(Y2)[-length(colnames(Y2))]
  rownames(A) <- colnames(Y2)[-length(colnames(Y2))]
}else{
  
  for(i in 1:dim(P)[1]){
    Y2 <- data.frame(t(Y))

    # infer with lasso
    fit <- glmnet(Y2, as.matrix(-P[i,]), alpha = 1)
    tmpA <- as.matrix(coef(fit, s = eps(1)))[-1,]
    A[,i] <- unname(tmpA)
    
    ## perm test
    tmpAS <- list()
    for (p in 1:permn){
      fitS <- glmnet(Y2, sample(as.matrix(-P[i,])), alpha = 1)
      tmpAS[[p]] <- as.matrix(coef(fitS, s = eps(1)))[-1,]
    }
    tmpAS2 <- do.call(rbind.data.frame, tmpAS)
    p <- c()
    for(k in 1:length(tmpA)){
      if(tmpA[k] <= 0){
        p[k] <- (length(which(tmpAS2[,k] <= as.numeric(tmpA[k])))/permn)+1/permn
      }else
      {
        p[k] <- (length(which(tmpAS2[,k] >= as.numeric(tmpA[k])))/permn)+1/permn
      }
    }
    A[match(names(tmpA)[which(p >= pval)],colnames(Y2)),i] <- 0
    progress(i,dim(P)[1], progress.bar = T)
    if(i == dim(P)[1]) message("Done!")
  }
  diag(A) <- 0
  A <- A %>% as.data.frame
  colnames(A) <- colnames(Y2)
  rownames(A) <- colnames(Y2)
}
return(A)
}