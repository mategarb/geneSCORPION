intersectNets <- function(Net1, Net2, diag=FALSE){
  
  if(diag == FALSE)
  {
    diag(Net1) <- 0
    diag(Net2) <- 0
  }
  Net1[Net1 < 0] <- -1
  Net1[Net1 > 0] <- 1
  
  Net2[Net2 < 0] <- -1
  Net2[Net2 > 0] <- 1
  
  Net3 <- (abs(Net1) + abs(Net2))/2
  Net3[Net3 < 1] <- 0
  Net4 <- zeros(dim(Net1)[1])
  Net4[Net3 == 1] <- Net1[Net3 == 1] 
  return(Net4)
}