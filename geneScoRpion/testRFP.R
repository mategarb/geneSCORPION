
library(R.matlab)
library(tidyverse)
library(pracma)
library(igraph)

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
  Anet <- readMat(paste0("/Users/mateuszgarbulowski/Desktop/datanets2/Anet_g50_S3_snr0001_rep2_",i,".mat"))
  Anet <- Anet$outA
  
  Ye <- read.csv2(file="/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/datasets/ymatrix_hepg2.csv", header=T, sep="\t",row.names = 1)
  Ye <- as.matrix(sapply(Ye, as.numeric))  
  #Y <- readMat(paste0("/Users/mateuszgarbulowski/Desktop/datanets2/y_g50_S3_snr0001_rep2_",i,".mat")) 
  #Y <- Y$outY
  P <- matrix(0, dim(Ye)[1], dim(Ye)[1])
  P <- diag(dim(Ye)[1])
  P <- cbind(P,P)
  
  
  out1 <- geneRanForest(Ye, P, 500, 1000)
  saveRDS(out1, "/Users/mateuszgarbulowski/Desktop/RFPerms_500tx1000p.rds")
  nlinks <- c()
  nlinkspa <- c()
  ps <- c(0, 0.001, 0.005, 0.01, 0.05, 0.1)
  S <- Spa <- c()
  for(i in 1:length(ps)){
  pvs <- out1$pValues
  pvs2 <- pvs <= ps[i]
  diag(pvs2) <- F
  nlinks[i] <- sum(pvs2)
  S[i] <- mean(colSums(pvs2 != 0))
  
  pvs <- out1$pValues
  pvs2 <- pAdjustNet(pvs,"fdr","col.wise") <= ps[i]
  diag(pvs2) <- F
  nlinkspa[i] <- sum(pvs2)
  Spa[i] <- mean(colSums(pvs2 != 0))
  }
  
  library(plotly)
  pvalues <- ps
  nlinks <- nlinkspa
  sparsity <- Spa

  plot_ly(x=pvalues, y=nlinks, z=sparsity, type="scatter3d", mode="markers", color=pvalues) %>%
     layout(scene = list(xaxis = list(title = 'P value'),
                            yaxis = list(title = 'number of links'),
                            zaxis = list(title = 'avg. sparsity per gene')))

  
  TFs <- read.table("/Users/mateuszgarbulowski/Desktop/TF_names_v_1.01.txt")
  TFs <- as.character(as.matrix(TFs)) #http://humantfs.ccbr.utoronto.ca/

  # calculate sign from LSCO
  nets_lsco <- infeRanea(Ye, P, method="lsco", zetas=0)
  nets_lsco2 <- nets_lsco[[1]]
  nets_lsco2[nets_lsco2 < 0] <- -1
  nets_lsco2[nets_lsco2 > 0] <- 1
  
  pvs <- pAdjustNet(out1$pValues,"fdr","col.wise") <= 0.05
  A1 <- nets_lsco2
  A1[which(!pvs)] <- 0
  
  ## net  
  hepg2 <- read.csv2(file="/Users/mateuszgarbulowski/Desktop/Postdoc SU_SciLifeLab/datasets/ymatrix_hepg2.csv", header=T, sep="\t",row.names = 1)
  colnames(A1) <- rownames(hepg2)
  rownames(A1) <- rownames(hepg2)
  nets_lsco2 <- as.data.frame(nets_lsco2)
  colnames(nets_lsco2) <- rownames(hepg2)
  rownames(nets_lsco2) <- rownames(hepg2)
  net2 <- graph_from_adjacency_matrix(
    t(A1),
    mode = c("directed"),
    weighted = T)
  
  net3 <- simplify(net2, remove.multiple = T, remove.loops = T) 
  
  # construct edges
  edges <- as_data_frame(net3, what = c("edges"))
  edges$color <- rep("dimgray",length(edges$from))
  
  aShape <- as.factor(edges$weight)
  levels(aShape) <- c("bar","arrow")
  edges$arrows.from.enabled <- TRUE
  edges$arrows.from.type <- as.character(aShape)
  
  #construct nodes
  nodes <- data.frame(id=unique(c(edges$from,edges$to)), shape='circle')
  nodes$label <- unique(c(edges$from,edges$to))
  nodes$color <- rep("darkturquoise",length(nodes$label))
  nodes$color[which(nodes$id %in% TFs)] <- "tomato"
  nodes$shape <- "dot"
  
  
  ## CYTOSCAPE
  inter <- ones(1,length(edges$weight))
  inter[edges$weight<0] <- -1
  inter <- as.factor(inter)
  levels(inter) <- c("inhibits","activates")
  inter2 <- as.character(inter)
  
  group <- rep("TG", length(nodes$id))
  group[which(nodes$id %in% TFs)] <- "TF"
  
  library(RCy3)
  style.name = "myStyle1"
  defaults <- list(NODE_SHAPE="ellipse",
                   NODE_SIZE=20,
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="W,E,c,0.00,0.00")
  nodeLabels <- mapVisualProperty('node label','id','p')
  nodeFills <- mapVisualProperty('node fill color','group','d',c("TG","TF","P"), c("slategray","orangered","cyan"))
  arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
  #edgeWidth <- mapVisualProperty('edge width','p')
  
  createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes))
  
  
  nodes2 <- data.frame(id=nodes$id,
                       #color=nodes$color, # categorical strings
                       group=group,
                       stringsAsFactors=FALSE)
  edges2 <- data.frame(source=edges$from,
                       target=edges$to,
                       interaction=inter2,  # optional
                       weight=edges$weight, # numeric
                       stringsAsFactors=FALSE)
  createNetworkFromDataFrames(nodes2,edges2, title="hepg2", collection="ENCODE")
  setVisualStyle(style.name)
  
  hubs <- hub_score(net3, scale=T)
  #sort(hubs$vector)
  
  tophubs <- names(sort(hubs$vector, decreasing = T)[1:10])
 
  
  edges3 <- edges2[which(edges2$source %in% tophubs),]
  nodes3 <- nodes2[match(unique(c(edges3$source,edges3$target)),nodes2$id),]
  createNetworkFromDataFrames(nodes3,edges3, title="hepg2", collection="ENCODE")
  setVisualStyle(style.name)
  
  alltfs <- nodes2$id[which(nodes2$group=="TF")]
  edges3 <- edges2[which(edges2$target %in% alltfs),]
  nodes3 <- nodes2[match(unique(c(edges3$source,edges3$target)),nodes2$id),]
  createNetworkFromDataFrames(nodes3,edges3, title="hepg2", collection="ENCODE")
  setVisualStyle(style.name)
  
  
  sq <- seq(0,1,0.01)
  outAUC <- c()
  for(i in 1:length(sq)){
  pVals3 <- out1$pValues <= sq[i]
  
  A2 <- matrix(0,ncol = dim(P)[1], nrow = dim(P)[1])
  A2[pVals3] <- 1
  outAUC[i] <- comp2Nets(Anet, A2, diag = T)$AUROC
  }
  #plot(sq,outAUC)
  
  nets_rf <- infeRanea(Y, P, method="rf", zetas = linspace(0,6,100))
  AUCsNoP <- lapply(nets_rf, function(x) comp2Nets(Anet,x,diag = F)$AUROC) %>% unlist

  #plot(linspace(0,6,100),AUCsNoP)
  
  
  par(mar=c(5,5,5,5)+0.1, las=1)
  
  plot.new()
  plot.window(xlim=range(sq), ylim=range(outAUC))
  points(sq, outAUC, col="red", pch=19)
  axis(1)
  axis(2, col.axis="red")
  box()
  
  plot.window(xlim=range(linspace(0,6,100)), ylim=range(AUCsNoP))
  points(linspace(0,6,100), AUCsNoP, col="limegreen", pch=19)
  axis(4, col.axis="limegreen")
  
  title("SNR = 0.001", adj=0)
  mtext("AUC (zeta linear)", side = 4, las=3, line=3, col="limegreen")
  mtext("AUC (p-value)", side = 2, las=3, line=3, col="red")
  