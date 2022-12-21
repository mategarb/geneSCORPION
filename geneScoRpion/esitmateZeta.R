library(R.matlab)
library(tidyverse)
library(pracma)

source("/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/R packages/geneRanea/comp2Nets.R")
source("/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/R packages/geneRanea/infeRanea.R")
source("/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/R packages/geneRanea/netwoRanea.R")
source("/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/R packages/geneRanea/spsity.R")

spsityLst <- function(NET){
  S1 <- S2 <- c()
  for(i in 1:length(NET)){
    NET2 <- NET[[i]]
    S1[i] <- spsity(NET2, "pergene")
    S2[i] <- spsity(NET2, "overall")
  }
  return(list(S1,S2))
}

AUCs1 <- AUCs2 <- AUCs3 <- AUCs4 <- AUCs5 <- AUCs6 <- list()
MCCs1 <- MCCs2 <- MCCs3 <- MCCs4 <- MCCs5 <- MCCs6 <- list()
F11 <- F12 <- F13 <- F14 <- F15 <- F16 <- list() 
  i=1
  Anet <- readMat(paste0("/Users/mateuszgarbulowski/Desktop/datanets2/Anet_g50_S3_snr001_rep2_",i,".mat"))
  Anet <- Anet$outA
  
  #Ye <- read.csv2(file="/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/datasets/ymatrix_hepg2.csv", header=T, sep="\t",row.names = 1)
  #Ye <- as.matrix(sapply(Ye, as.numeric))  
  Y <- readMat(paste0("/Users/mateuszgarbulowski/Desktop/datanets2/y_g50_S3_snr001_rep2_",i,".mat")) 
  Y <- Y$outY
  P <- matrix(0, dim(Y)[1], dim(Y)[1])
  P <- diag(dim(Y)[1])
  P <- cbind(P,P)
  
  ## lsco
  # starting point - draw ranges
  deSpa <- 10
  
  minT1 <- runif(n=1, min = 0, max = 1)
  maxT1 <- runif(n=1, min = minT1, max = 1)
  thr1 <- runif(n=1, min = minT1, max = maxT1)
  
  minT2 <- runif(n=1, min = 0, max = 1)
  maxT2 <- runif(n=1, min = minT2, max = 1)
  thr2 <- runif(n=1, min = minT2, max = maxT2)
  
  keepMe <- keepMe2 <- c()
  for(i in 1:100){
  nets1 <- infeRanea(Y, P, method="lsco", zetas=c(thr1, thr2))
  wm <- which.min(abs(spsityLst(nets1)[[1]]-deSpa))
  
  if(wm == 1){
    minT2 <- runif(n=1, min = 0, max = 1)
    maxT2 <- runif(n=1, min = minT2, max = 1)
    thr2 <- runif(n=1, min = minT2, max = maxT2)
  }else{
    minT1 <- runif(n=1, min = 0, max = 1)
    maxT1 <- runif(n=1, min = minT1, max = 1)
    thr1 <- runif(n=1, min = minT1, max = maxT1)
  }
  keepMe[i] <- spsityLst(nets1)[[1]][wm]
  keepMe2[i] <-  c(thr1,thr2)[wm]
  #nets2 <- infeRanea(Y, P, method="lsco", zetas=c(thr1, thr2))
  #wm2 <- which.min(abs(spsityLst(nets2)[[1]]-deSpa))
  
  
  if (spsityLst(nets1)[[1]][wm] == deSpa) {
    break
  }
  }
  
  plot(keepMe2, keepMe, type="l")
  
  nets2 <- infeRanea(Y, P, method="lsco", zetas=c(thr1, thr2))
  wm2 <- which.min(abs(spsityLst(nets2)[[1]]-deSpa))
  
  c(thr1, thr2)[wm2]
  
  ### AUROC
  
  deAUC <- 1
  
  minT1 <- runif(n=1, min = 0, max = 1)
  maxT1 <- runif(n=1, min = minT1, max = 1)
  thr1 <- runif(n=1, min = minT1, max = maxT1)
  
  minT2 <- runif(n=1, min = 0, max = 1)
  maxT2 <- runif(n=1, min = minT2, max = 1)
  thr2 <- runif(n=1, min = minT2, max = maxT2)
  
  keepMe <- keepMe2 <- c()
  for(i in 1:100){
    nets1 <- infeRanea(Y, P, method="lsco", zetas=c(thr1, thr2))
    wm <- which.min(abs(unlist(lapply(nets1, function(x) comp2Nets(Anet,x,diag = F)$AUROC))-deAUC))
    
    if(wm == 1){
      minT2 <- runif(n=1, min = 0, max = 1)
      maxT2 <- runif(n=1, min = minT2, max = 1)
      thr2 <- runif(n=1, min = minT2, max = maxT2)
    }else{
      minT1 <- runif(n=1, min = 0, max = 1)
      maxT1 <- runif(n=1, min = minT1, max = 1)
      thr1 <- runif(n=1, min = minT1, max = maxT1)
    }
    keepMe[i] <- unlist(lapply(nets1, function(x) comp2Nets(Anet,x,diag = F)$AUROC))[wm]
    keepMe2[i] <-  c(thr1,thr2)[wm]
    #nets2 <- infeRanea(Y, P, method="lsco", zetas=c(thr1, thr2))
    #wm2 <- which.min(abs(spsityLst(nets2)[[1]]-deSpa))
    
    
    if (unlist(lapply(nets1, function(x) comp2Nets(Anet,x,diag = F)$AUROC))[wm] == deAUC) {
      break
    }
  }
  
  plot(keepMe2, keepMe, type="l")
  c(thr1, thr2)[wm]
  