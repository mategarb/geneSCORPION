pAdjustNet <- function(pValues, method, option="all"){
  
  if(option == "all"){
    adjpVals <- matrix(p.adjust(as.vector(as.matrix(pValues)), method=method), ncol=dim(pValues)[1])
  }  
  
  if(option == "row.wise"){
    adjpVals <- apply(pValues, 1, function(x) p.adjust(x, method))
  }  
  
  if(option == "col.wise"){
    adjpVals <- apply(pValues, 2, function(x) p.adjust(x, method))
  }
  
  return(adjpVals)
}

