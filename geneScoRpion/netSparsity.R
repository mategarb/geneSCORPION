netSparsity <- function(Net, type="overall"){
  if(type=="overall"){
    S <- length(which(Net!=0))/(dim(Net)[1]*dim(Net)[2])
  }
  if(type=="pergene"){
    S <- mean(colSums(Net != 0))
  }
  return(S)
}