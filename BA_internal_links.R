# The implementation of the Barabasi-Albert Model allowing for internal links. 
# The initial number of nodes is fixed at 2.


BA_internal <- function(N, m = 2, n = 0){
  
  m_0 <- 2 #initial nodes
  t <- N - m_0 #time steps
  nodes <- c(1,2)
  edges <- tibble(from = 1, to = 2)
  for (s in 1:t){
    new_node <- max(nodes) + 1 #new node
    nodes_d <- NULL
    for (j in nodes){
      nodes_d[j] <- sum(edges$from == j) + sum(edges$to == j)
      - sum((edges$from == j) & (edges$to == j))
    }
    
    #m new links are created
    if (length(nodes) >= m){
      a <- sample(nodes, m, replace = F, prob = nodes_d)
      b <- rep(new_node, m)
    }else{
      a <- sample(nodes, length(nodes), replace = F, prob = nodes_d) #PA
      b <- rep(new_node, length(nodes))
    }
    edges <- rbind(edges, tibble(from = b, to = a))
    nodes <- c(nodes, new_node)
    nodes_d <- NULL
    for (j in nodes){
      nodes_d[j] <- sum(edges$from == j) + sum(edges$to == j)
      - sum((edges$from == j) & (edges$to == j))
    }
    
    #After the new node has been connected, if the number of internal
    #links n is greater than 0, n new links are generated
    if (n > 0){
      nodes_couple <- NULL
      nodes_couple_d <- NULL
      for (i in nodes){
        for (j in i:max(nodes)){
          nodes_couple <- c(nodes_couple, list(c(i, j)))
          nodes_couple_d <- c(nodes_couple_d, nodes_d[i]*nodes_d[j]) #double PA
        }
      }
      for (int in 1:n){
        a <- sample(nodes_couple, 1, prob = nodes_couple_d)
        edges <- rbind(edges, c(a[[1]][1], a[[1]][2]))
      }
    }
  }
  #Returning the network
  network <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
  network
}