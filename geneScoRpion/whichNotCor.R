whichNotCor <- function(Y){
  Y <- Y %>% as.data.frame
  vec <- c()
  vec2 <- c()
  for(i in 1:dim(Y)[1]){
    
    cr <- cor(Y[,c(i,dim(Y)[1]+i)])
    pcr <- cor.test(Y[,c(i,dim(Y)[1]+i)][,1],Y[,c(i,dim(Y)[1]+i)][,2])
    vec[i] <- cr[1,2]
    vec2[i] <- pcr$p.value
  }
  remCols <- colnames(Y)[which(p.adjust(vec2,method = "fdr")>=0.05)]
  colInd <- which(grepl(paste0(remCols,collapse = "|"),colnames(Y)))
  rowInd <- which(grepl(paste0(remCols,collapse = "|"),rownames(Y)))
  return(list(colInd,rowInd))
}