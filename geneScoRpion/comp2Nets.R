comp2Nets <- function(Net1, Net2, diag=TRUE) {
#nets as numeric matrices of the same dimensions
tab0 <- table(c(0,1),c(1,0))
tab0[1,2] <- 0
tab0[2,1] <- 0

if(diag==FALSE)
{
  diag(Net1) <- 0
  diag(Net2) <- 0
}
rownames(Net1) <- NULL
colnames(Net1) <- NULL
rownames(Net2) <- NULL
colnames(Net2) <- NULL
Net1 <- as.matrix(Net1)
Net2 <- as.matrix(Net2)
for(i in 1:dim(Net1)[1]){
  Net1[i,which(Net1[i,] != 0)] <- 1
  Net2[i,which(Net2[i,] != 0)] <- 1
  v1 <- factor(Net1[i,])
  v2 <- factor(Net2[i,])
  levels(v1) <- c(0,1)
  levels(v2) <- c(0,1)
  tab0 <- tab0+table(v1,v2)
}
n1 <- as.numeric(Net1)
n2 <- as.numeric(Net2)

auroc <- suppressMessages(auc(n1, n2))

#return confusion matrix
colnames(tab0) <- c("Negative","Positive")
rownames(tab0) <- c("Negative","Positive")

#c("TN","FN")
#c("FP","TP")
TN <- tab0[1,1]
TP <- tab0[2,2]
FP <- tab0[2,1]
FN <- tab0[1,2]
ACC <- (TP+TN)/(TP+FP+FN+TN)
rec <- TP/(TP+FN)
prec <- TP/(TP+FP)
F1 <- 2*(prec*rec)/(prec+rec)
MCC <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
S1 <- netSparsity(Net1,"pergene")
S2 <- netSparsity(Net2,"pergene")
nlinks <- length(which(Net2!=0))
return(list('TN'=TN,'TP'=TP,'FP'=FP,'FN'=FN, 'ACC'=ACC, 'F1'=F1, 'MCC'=MCC,'S1'=S1, 'S2'=S2,'nlinks'=nlinks, 'AUROC'=auroc))
}