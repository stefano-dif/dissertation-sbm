# The implementation of the DC-SBM generation algorithm, undirected and directed.
# It allows to set a parameter lambda to interpolate between a random graph and
# any specified community structure. If no community structure is specified, an assortative
# one is chosen. The elements of the group vector are expected to range from 1 to the number of communities.


# DCSBM - Undirected

DCSBM_und <- function(n, g, d, omega_planted = "canonical", lambda = 0.5){
  
  m <- sum(d)/2 #number of edges
  q <- max(g) #number of communities
  k <- NULL
  for (i in 1:q){
    k[i] <- sum(d[g == i]) #kappas
  }
  theta <- d/(k[g]) #theta vector
  #Matrix Omega:
  omega_random <- matrix(nrow = q, ncol = q)
  for (i in 1:q){
    for (j in 1:q){
      omega_random[i, j] <- (k[i]*k[j])/(2*m)
    }
  }
  if (all(omega_planted == "canonical")){
    omega_planted <- matrix(nrow = q, ncol = q)
    for (i in 1:q){
      for (j in 1:q){
        if (i == j){
          omega_planted[i, j] <- k[i]
        }else{
          omega_planted[i, j] <- 0
        }
      }
    }
  }
  omega <- lambda*omega_planted + (1 - lambda)*omega_random
  
  #Step 1: drawing the total number of edges between two groups
  #from a Poisson distribution
  ms <- matrix(nrow = q, ncol = q)
  for (i in 1:q){
    for (j in i:q){
      if (i == j){
        ms[i, j] <- rpois(1, omega[i, j]/2)
      }else{
        ms[i, j] <- rpois(1, omega[i, j])
      }
    }
  }
  
  #Step 2: connecting each edge to dyads according to theta
  nodes <- 1:n
  edges <- NULL
  for (g1 in 1:q){
    for (g2 in g1:q){
      if (ms[g1, g2] > 0){
        for (e in 1:ms[g1, g2]){
          a <- sample(nodes[g == g1], 1, prob = theta[g == g1])
          b <- sample(nodes[g == g2], 1, prob = theta[g == g2])
          temp <- c(a, b)
          edges <- rbind(edges, temp)
        }
      }
    }
  }
  
  #Returning the network
  edges <- tibble(from = edges[, 1], to = edges[, 2])
  rownames(edges) <- NULL
  nodes <- tibble(id = 1:n, group = g)
  network <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
}


# DCSBM - Directed

DCSBM_D <- function(n, g, d = "ab", omega_planted = "canonical", lambda = 0.5){
  if (all(d == "ab")){
    base <- sample_pa(n, power = 1)
    d <- degree(base)
    m <- length(E(base))
  }else{
    m <- sum(d)/2
  }
  q <- max(g)
  k <- NULL
  for (i in 1:q){
    k[i] <- sum(d[g == i])
  }
  theta <- d/(k[g])
  
  omega_random <- matrix(nrow = q, ncol = q)
  for (i in 1:q){
    for (j in 1:q){
      omega_random[i, j] <- (k[i]*k[j])/(2*m)
    }
  }
  
  if (all(omega_planted == "canonical")){
    omega_planted <- matrix(nrow = q, ncol = q)
    for (i in 1:q){
      for (j in 1:q){
        if (i == j){
          omega_planted[i, j] <- k[i]
        }else{
          omega_planted[i, j] <- 0
        }
      }
    }
  }
  
  omega <- lambda*omega_planted + (1 - lambda)*omega_random
  
  
  
  #Step 1:
  
  ms <- matrix(nrow = q, ncol = q)
  
  for (i in 1:q){
    for (j in 1:q){
      if (i == j){
        ms[i, j] <- rpois(1, omega[i, j]/2)
      }else{
        ms[i, j] <- rpois(1, omega[i, j])
      }
    }
  }
  
  
  #Step 2:
  
  nodes <- 1:n
  edges <- NULL
  
  for (g1 in 1:q){
    for (g2 in 1:q){
      if (ms[g1, g2] > 0){
        for (e in 1:ms[g1, g2]){
          a <- sample(nodes[g == g1], 1, prob = theta[g == g1])
          b <- sample(nodes[g == g2], 1, prob = theta[g == g2])
          temp <- c(a, b)
          edges <- rbind(edges, temp)
        }
      }
    }
  }
  
  edges <- tibble(from = edges[, 1], to = edges[, 2])
  rownames(edges) <- NULL
  
  nodes <- tibble(id = 1:n, group = g)
  
  network <- graph_from_data_frame(d = edges, vertices = nodes, directed = T)
  
  
}