pertRank <- function(Y, P){
  rY <- apply(Y,1,rank)/dim(Y)[2]
  rY2 <- rY[P == 1]
  return(rY2)
}