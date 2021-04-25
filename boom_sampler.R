### Sampler

boom_sampler <- function(A,
                         G,
                         y,
                         S){
  #############################################################################
  # Set-Up
  #############################################################################
  # Problem Dimensions
  n <- length(y)
  V <- dim(G)[2]
  P <- dim(G)[3]
  
  # Data Pre-processing
  X      <- G
  dim(X) <- c(n, P * V)
  y      <- y
  
  # Variable Initialization
  xi <- rep(1, P)
  v  <- matrix(data = 1, nrow = V, ncol = P)
  t2 <- rep(1, P)
  l2 <- matrix(data = 1, nrow = V, ncol = P)
  s2 <- 1
  g  <- rep(1, P)
  e2 <- rep(1, P)
  
  # Sample Variables
  sam_B     <- array(data = NA, dim = c(S, V, P))
  sam_l2    <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2 <- array(data = NA, dim = c(S, V, P))
  sam_v     <- array(data = NA, dim = c(S, V, P))
  sam_g     <- matrix(data = NA, nrow = S, ncol = P)
  sam_pg    <- matrix(data = NA, nrow = S, ncol = P)
  sam_e2    <- matrix(data = NA, nrow = S, ncol = P)
  sam_t2    <- matrix(data = NA, nrow = S, ncol = P)
  sam_xi    <- matrix(data = NA, nrow = S, ncol = P)
  sam_s2    <- numeric(length = S)
  
  #############################################################################
  # Sampling
  #############################################################################
  # Progress Bar
  pb <- txtProgressBar(min     = 0,
                       max     = 1,
                       initial = 0,
                       style   = 3,
                       width   = 72)
  
  # Sampling
  for(s in 1:S){
    ###########################################################################
    # Samples B Structure
    ###########################################################################
    # Samples the Horseshoe Structure
    out <- fast_spike_horseshoe(y  = y,
                                X  = X,
                                xi = xi,
                                v  = v,
                                l2 = l2,
                                t2 = t2,
                                s2 = s2,
                                g  = g,
                                e2 = e2)
    # Assigns the values
    B  <- out$b
    s2 <- out$s2
    l2 <- out$l2
    t2 <- out$t2
    v  <- out$v
    xi <- out$xi
    e2 <- out$e2
    g  <- out$g
    pg <- out$pg
    
    # Saves Variables
    sam_B[s,,]     <- B
    sam_l2[s,,]    <- l2
    sam_l2_t2[s,,] <- t(t(l2) * t2)
    sam_v[s,,]     <- v
    sam_t2[s,]     <- t2
    sam_xi[s,]     <- xi
    sam_s2[s]      <- s2
    sam_e2[s,]     <- e2
    sam_g[s,]      <- g
    sam_pg[s,]     <- pg
    
    ###########################################################################
    # Progress Bar Update
    ###########################################################################
    setTxtProgressBar(pb    = pb,
                      value = s / S)
    
  }
  
  # Returns
  return(list(B     = sam_B,
              s2    = sam_s2,
              l2    = sam_l2,
              t2    = sam_t2,
              v     = sam_v,
              xi    = sam_xi,
              g     = sam_g,
              e2    = sam_e2,
              l2_t2 = sam_l2_t2,
              pg    = sam_pg))
}