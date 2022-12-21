geneRanForest <- function(Y, P, treesNum=100, perms=10){
    
    # load libraries
    library(ranger)
  
    # initialize progress bar
    pb <- txtProgressBar(min = 0, max = dim(P)[1], style = 3, width = 50, char = "=")
    
    # create empty matrix to store importance from trees
    A <- matrix(0,ncol = dim(P)[1], nrow = dim(P)[1])
    
    # create an empty list of matrices for permutation procedure
    Aperm <- list()
    for(g in 1:perms){
      Aperm[[g]] <- matrix(0,ncol = dim(P)[1], nrow = dim(P)[1])
    }
    
    # run loop for all the genes
    for(i in 1:dim(P)[1]){
      
      # infer a network for non-permuted data
      data <- data.frame(t(Y), decision=unname(t(P[i,])))
      rf <- ranger(decision~., data, importance = "impurity", probability = T,
                   classification = T, num.trees = treesNum, splitrule = "extratrees", replace = F)
      A[i,] <- unname(rf$variable.importance)
      
      # run permutations by shuffling decision, i.e. perturbed vs non-perturbed
      for(j in 1:perms){
        dataPerm <- data.frame(t(Y), decision=unname(sample(as.matrix(P[i,]))))
        rfPerm <- ranger(decision~., dataPerm, importance = "impurity", probability = T,
                     classification = T, num.trees = treesNum, splitrule = "extratrees", replace = F)
        Aperm[[j]][i,] <- unname(rfPerm$variable.importance)
      }
      setTxtProgressBar(pb, i) # update progress bar
    }
    
    # create an empty matrix for P values
    pVals <- matrix(0,ncol = dim(P)[1], nrow = dim(P)[1])
    ## calculate p-values by comparing original value to random distribution and divide by the number of permutations
    for(i in 1:dim(P)[1]){
      for(j in 1:dim(P)[1]){
        pVals[i,j] <- sum(as.numeric(lapply(Aperm, function(x) x[i,j])) > A[i,j])/perms
      }
    }
    
    # output a list that contains infered network and P values 
    return(list(fullNet=t(A), pValues=t(pVals)))

}
