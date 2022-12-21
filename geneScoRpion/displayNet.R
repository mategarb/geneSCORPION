displayNet <- function(A){
net2 <- graph_from_adjacency_matrix(
  as.matrix(A),
  mode = c("directed"),
  weighted = T)

net3 <- simplify(net2, remove.multiple = T, remove.loops = T) 

# construct edges
edges <- as_data_frame(net3, what = c("edges"))
edges$dashes <- rep(FALSE, length(edges$from))
#edges$dashes[which(edges$weight<0)]<- TRUE
#edges$weight <- c()
edges$arrows <- rep("from", length(edges$from))
edges$color <- rep("dimgray",length(edges$from))

#construct nodes
nodes <- data.frame(id=unique(c(edges$from,edges$to)), shape='circle')
nodes$label <- unique(c(edges$from,edges$to))
nodes$color <- rep("darkturquoise",length(nodes$label))

visNetwork(nodes, edges) %>% visLayout(randomSeed = 1)
}