# The implementation of the DCSBM-based and SBM-based Karrer and Newman algorithm, together with their helper functions.
# It can be generalized to all methods seen in Section 2 of the thesis by changing the function computing the likelihood.


#DCSBM-based


# mrs: number of edges between two groups r and s

mrs <- function(network, r, s){
  #adjacency matrix of the network
  adj <- as_adj(network, sparse = F)[V(network)$community == r,
                                     V(network)$community == s]
  temp <- sum(adj) #number of edges between two distinct groups
  
  #In the case of r = s:
  if (r == s){
    loops <- sum(diag(adj)) #number of self-loops
    temp <- (temp - loops)/2 + loops
    temp <- temp * 2
    #notice that as_adj counts the number of edges even on the main diagonal
  }
  result <- temp
  result
}


#kappa: sum of degrees of nodes in a group r

kappa <- function(network, r){
  sum(degree(network)[V(network)[V(network)$community == r]])
}


#Log-Likelihood (DC-SBM)

likelihood <- function(network){
  l <- 0
  comms <- unique(V(network)$community)
  for (i in comms){
    for (j in comms){
      if(mrs(network, i, j) == 0){
        l <- l + 0
      }else{
        l <- l +
          mrs(network, i, j)*log(mrs(network, i, j)/(kappa(network, i)*kappa(network, j)))
      }
    }
  }
  l
}



#DCSBM Node-moving algorithm

CD_DCSBM <- function(network, q, c = "random") {
  
  iter <- 1 #Iteration number
  if(all(c == "random")){
    c <- sample(1:q, size = length(V(network)), replace = T)
  } #Set the initial partition
  V(network)$community <- c
  cat("Initial Likelihood:", likelihood(network), "\n")
  
  #Algorithm iteration
  repeat{
    cat("Iteration:", iter, "\n")
    #Benchmark (likelihood)
    if(iter == 1){
      benchmark_l <- likelihood(network)
    }else{
      benchmark_l <- likes[best_state]
    }
    moved <- length(V(network)) + 1 #to store moved nodes
    likes <- likelihood(network) #to store the likelihood of each state
    inspection_table <- NULL #to implement the best state
    
    #Repeating the step for each node in the network
    for(step in 1:length(V(network))){
      cat("Step:", step, " ")
      options <- NULL #vector of the likelihood deltas
      register <- NULL #store moves
      
      #Find the best possible move among nodes which have not been moved yet
      for (i in V(network)[-moved]){
        ir <- V(network)$community[i]
        trials <- unique(V(network)$community)[unique(V(network)$community) != ir]
        for (s in trials){
          r <- V(network)$community[i]
          V(network)$community[i] <- s
          new_like <- likelihood(network)
          V(network)$community[i] <- r
          options <- c(options, new_like - benchmark_l)
          reg <- data.frame(node = i, to = s)
          register <- rbind(register, reg)
        }
      }
      if (length(which(options == max(options))) > 1){
        mx <- sample(which(options == max(options)), 1)
      }else{
        mx <- which(options == max(options))
      }
      moved_id <- register[mx, ]$node
      new_comm <- register[mx, ]$to
      action <- tibble(id = moved_id, from = V(network)$community[moved_id])
      #Implement the best move
      V(network)$community[moved_id] <- new_comm
      moved <- c(moved, moved_id)
      likes <- c(likes, likelihood(network))
      inspection_table <- rbind(inspection_table, action)
    }
    if (length(which(likes == max(likes))) > 1){
      best_state <- sample(which(likes == max(likes)), 1)
    }else{
      best_state <- which(likes == max(likes))
    }
    #Implement the state where the system reached maximum likelihood
    if (best_state <= nrow(inspection_table)){
      temp <- inspection_table[best_state:nrow(inspection_table), ]
      V(network)[temp$id]$community <- temp$from
    }
    #Condition terminating the algorithm
    if(likes[best_state] <= benchmark_l){
      #Algorithm ends
      cat("Stopped.")
      break
    }
    iter <- iter + 1
    cat("Likelihood:", likelihood(network), "\n")
  }
  #Returning the network with attribute 'community' containing the recovered partition
  return(network)
}


#----


# SBM-based


n <- function(network, r){
  sum(V(network)$community == r)
}



likelihood_sbm <- function(network){
  l <- 0
  comms <- unique(V(network)$community)
  for (i in comms){
    for (j in comms){
      if(mrs(network, i, j) == 0){
        l <- l + 0
      }else{
        l <- l + 
          mrs(network, i, j)*log(mrs(network, i, j)/(n(network, i)*n(network, j)))
        
      }
    }
  }
  l
}



CD_SBM <- function(network, q, c = "random") {
  
  iter <- 1
  
  if(all(c == "random")){
    c <- sample(1:q, size = length(V(network)), replace = T)
  }
  
  V(network)$community <- c
  cat("Initial Likelihood:", likelihood_sbm(network), "\n")
  
  repeat{
    cat("Iteration:", iter, "\n")
    
    if(iter == 1){
      benchmark_l <- likelihood_sbm(network)
    }else{
      benchmark_l <- likes[best_state]
    }
    
    moved <- length(V(network)) + 1 #to avoid the null case
    likes <- likelihood_sbm(network)
    inspection_table <- NULL
    
    for(step in 1:length(V(network))){
      cat("Step:", step, " ")
      
      options <- NULL 
      register <- NULL 
      
      for (i in V(network)[-moved]){
        ir <- V(network)$community[i]
        trials <- unique(V(network)$community)[unique(V(network)$community) != ir]
        for (s in trials){
          x <- V(network)$community[i]
          V(network)$community[i] <- s
          l2 <- likelihood_sbm(network)
          V(network)$community[i] <- x
          options <- c(options, l2 - benchmark_l)
          reg <- data.frame(node = i, to = s)
          register <- rbind(register, reg)
        }
      }
      
      if (length(which(options == max(options))) > 1){
        mx <- sample(which(options == max(options)), 1)
      }else{
        mx <- which(options == max(options))
      }
      
      moved_id <- register[mx, ]$node
      new_comm <- register[mx, ]$to
      action <- tibble(id = moved_id, from = V(network)$community[moved_id])
      
      V(network)$community[moved_id] <- new_comm
      moved <- c(moved, moved_id)
      likes <- c(likes, likelihood_sbm(network))
      inspection_table <- rbind(inspection_table, action)
    }
    
    if (length(which(likes == max(likes))) > 1){
      best_state <- sample(which(likes == max(likes)), 1)
    }else{
      best_state <- which(likes == max(likes)) 
    }
    
    if (best_state <= nrow(inspection_table)){
      temp <- inspection_table[best_state:nrow(inspection_table), ]
      V(network)[temp$id]$community <- temp$from
    }
    
    if(likes[best_state] <= benchmark_l){
      #Algorithm ends
      cat("Stopped.")
      break
    }
    
    iter <- iter + 1
    cat("Likelihood:", likelihood_sbm(network), "\n")
  }
  
  return(network)
}

#---------

# VARIANTS



log_ <- function(x){
  if (x != 0){
    log(x)
  }else{
    0
  }
}


fi <- function(k, f){
  max(f, 1/k)
}


H <- function(network, i, f){
  a <- kit(network, i, V(network)$community[i])
  b <- fi(degree(network)[i], f)
  -a*log_(b) -(1-a)*log_(1-b)
}


likelihood_ASSOR <- function(network, f){
  l1 <- 0
  comms <- unique(V(network)$community)
  for (i in comms){
    for (j in comms){
      if(mrs(network, i, j) == 0){
        l1 <- l1 + 0
      }else{
        l1 <- l1 + 
          mrs(network, i, j)*log(mrs(network, i, j)/(kappa(network, i)*kappa(network, j)))
        
      }
    }
  }
  l2 <- 0
  for (i in V(network)){
    l2 <- l2 + degree(network)[i]*H(network, i, f)
  }
  
  l3 <- 0
  for (i in V(network)){
    l3 <- l3 + log_(1 - fi(degree(network)[i], f))
  }
  
  l1 - 2*l2 + 2*l3
}


CD_DCSBM_reg <- function(network, q, f, c = "random") {
  
  iter <- 1
  
  if(all(c == "random")){
    c <- sample(1:q, size = length(V(network)), replace = T)
  }
  
  V(network)$community <- c
  cat("Initial Likelihood:", likelihood_ASSOR(network, f), "\n")
  
  repeat{
    cat("Iteration:", iter, "\n")
    
    if(iter == 1){
      benchmark_l <- likelihood_ASSOR(network, f)
    }else{
      benchmark_l <- likes[best_state]
    }
    
    moved <- length(V(network)) + 1 #to avoid the null case
    likes <- likelihood_ASSOR(network, f)
    inspection_table <- NULL
    
    for(step in 1:length(V(network))){
      cat("Step:", step, " ")
      
      options <- NULL 
      register <- NULL 
      
      for (i in V(network)[-moved]){
        ir <- V(network)$community[i]
        trials <- unique(V(network)$community)[unique(V(network)$community) != ir]
        for (s in trials){
          x <- V(network)$community[i]
          V(network)$community[i] <- s
          l2 <- likelihood_ASSOR(network, f)
          V(network)$community[i] <- x
          options <- c(options, l2 - benchmark_l)
          reg <- data.frame(node = i, to = s)
          register <- rbind(register, reg)
        }
      }
      
      if (length(which(options == max(options))) > 1){
        mx <- sample(which(options == max(options)), 1)
      }else{
        mx <- which(options == max(options))
      }
      
      moved_id <- register[mx, ]$node
      new_comm <- register[mx, ]$to
      action <- tibble(id = moved_id, from = V(network)$community[moved_id])
      
      V(network)$community[moved_id] <- new_comm
      moved <- c(moved, moved_id)
      likes <- c(likes, likelihood_ASSOR(network, f))
      inspection_table <- rbind(inspection_table, action)
    }
    #cat("\n")
    cat("likes_pre:", likes, "\n")
    
    likes[is.na(likes)] <- min(likes[!is.na(likes)]) - 1000 #correction
    
    if (length(which(likes == max(likes))) > 1){
      best_state <- sample(which(likes == max(likes)), 1)
    }else{
      best_state <- which(likes == max(likes)) 
    }
    
    if (best_state <= nrow(inspection_table)){
      temp <- inspection_table[best_state:nrow(inspection_table), ]
      V(network)[temp$id]$community <- temp$from
    }
    
    if(likes[best_state] <= benchmark_l){
      #Algorithm ends
      cat("Stopped.")
      break
    }
    
    iter <- iter + 1
    cat("Likelihood:", likelihood_ASSOR(network, f), "\n")
  }
  
  return(network)
}





kappa_out <- function(network, r){
  sum(degree(network, mode = "out")[V(network)[V(network)$community == r]])
}

kappa_in <- function(network, r){
  sum(degree(network, mode = "in")[V(network)[V(network)$community == r]])
}


likelihood_dir <- function(network){
  l <- 0
  comms <- unique(V(network)$community)
  for (i in comms){
    for (j in comms){
      if(mrs(network, i, j) == 0){
        l <- l + 0
      }else{
        l <- l + 
          mrs(network, i, j)*log(mrs(network, i, j)/(kappa_out(network, i)*kappa_in(network, j)))
        
      }
    }
  }
  l
}




CD_DCSBM_dir <- function(network, q, c = "random") {
  
  iter <- 1
  
  if(all(c == "random")){
    c <- sample(1:q, size = length(V(network)), replace = T)
  }
  
  V(network)$community <- c
  cat("Initial Likelihood:", likelihood_dir(network), "\n")
  
  repeat{
    cat("Iteration:", iter, "\n")
    
    if(iter == 1){
      benchmark_l <- likelihood_dir(network)
    }else{
      benchmark_l <- likes[best_state]
    }
    
    moved <- length(V(network)) + 1 #to avoid the null case
    likes <- likelihood_dir(network)
    inspection_table <- NULL
    
    for(step in 1:length(V(network))){
      cat("Step:", step, " ")
      
      options <- NULL 
      register <- NULL 
      
      for (i in V(network)[-moved]){
        ir <- V(network)$community[i]
        trials <- unique(V(network)$community)[unique(V(network)$community) != ir]
        for (s in trials){
          x <- V(network)$community[i]
          V(network)$community[i] <- s
          l2 <- likelihood_dir(network)
          V(network)$community[i] <- x
          options <- c(options, l2 - benchmark_l)
          reg <- data.frame(node = i, to = s)
          register <- rbind(register, reg)
        }
      }
      
      if (length(which(options == max(options))) > 1){
        mx <- sample(which(options == max(options)), 1)
      }else{
        mx <- which(options == max(options))
      }
      
      moved_id <- register[mx, ]$node
      new_comm <- register[mx, ]$to
      action <- tibble(id = moved_id, from = V(network)$community[moved_id])
      
      V(network)$community[moved_id] <- new_comm
      moved <- c(moved, moved_id)
      likes <- c(likes, likelihood_dir(network))
      inspection_table <- rbind(inspection_table, action)
    }
    #cat("\n")
    cat("likes_pre:", likes, "\n")
    
    likes[is.na(likes)] <- min(likes[!is.na(likes)]) - 1000 #correction
    
    if (length(which(likes == max(likes))) > 1){
      best_state <- sample(which(likes == max(likes)), 1)
    }else{
      best_state <- which(likes == max(likes)) 
    }
    
    if (best_state <= nrow(inspection_table)){
      temp <- inspection_table[best_state:nrow(inspection_table), ]
      V(network)[temp$id]$community <- temp$from
    }
    
    if(likes[best_state] <= benchmark_l){
      #Algorithm ends
      cat("Stopped.")
      break
    }
    
    iter <- iter + 1
    cat("Likelihood:", likelihood_dir(network), "\n")
  }
  
  return(network)
}



