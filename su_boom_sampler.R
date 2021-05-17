### Sampler

su_boom_sampler <- function(A,
                            G,
                            y,
                            S,
                            c2,
                            m2,
                            e2){
  #############################################################################
  # Set-Up
  #############################################################################
  # Problem Dimensions
  n <- dim(G)[1]
  V <- dim(G)[2]
  P <- dim(G)[3]
  
  # Data Pre-processing
  X1      <- G
  dim(X1) <- c(n, P * V)
  X2      <- A
  dim(X2) <- c(n, P * P)
  temp    <- matrix(data = NA, nrow = P, ncol = P)
  X2      <- X2[, c(lower.tri(temp))]
  X       <- cbind(X1, X2)
  
  # Variable Initialization
  l2       <- matrix(data = 1, nrow = V, ncol = P)
  v        <- matrix(data = 1, nrow = V, ncol = P)
  t2       <- 1
  xi       <- 1
  s2       <- 1
  g        <- c(rep(1, QQ), rep(0, P - QQ))
  #g        <- rep(0, P)
  
  # Sample Variables
  sam_B        <- array(data = NA, dim = c(S, V, P))
  sam_l2       <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2_s2 <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2    <- array(data = NA, dim = c(S, V, P))
  sam_e2_s2    <- array(data = NA, dim = c(S, V, P))
  sam_v        <- array(data = NA, dim = c(S, V, P))
  sam_t2       <- numeric(length = S)
  sam_xi       <- numeric(length = S)
  sam_Theta    <- array(data = NA, dim = c(S, P, P))
  sam_s2       <- numeric(length = S)
  sam_g        <- matrix(data = NA, nrow = S, ncol = P)
  sam_pg       <- matrix(data = NA, nrow = S, ncol = P)
  sam_D        <- matrix(data = NA, nrow = S, ncol = P * V + (P - 1) * P / 2)
  
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
    out <- su_boom_iteration(y  = y,
                             X  = X,
                             l2 = l2,
                             v  = v,
                             t2 = t2,
                             xi = xi,
                             c2 = c2,
                             m2 = m2,
                             s2 = s2,
                             g  = g,
                             e2 = e2)
    # Assigns and Distributes the values
    # B Structure
    B     <- out$B
    l2    <- out$l2
    v     <- out$v
    t2    <- out$t2
    xi    <- out$xi
    # Theta Structure
    Theta <- out$Theta
    # Error
    s2 <- out$s2
    # Sparsity Structure
    g  <- out$g
    pg <- out$pg
    # Coefficient Priors
    D  <- out$D
    
    # Saves Variables
    sam_B[s,,]        <- B
    sam_l2[s,,]       <- l2
    sam_l2_t2_s2[s,,] <- l2 * t2 * s2
    sam_l2_t2[s,,]    <- l2 * t2
    sam_e2_s2[s,,]    <- e2 * s2
    sam_v[s,,]        <- v
    sam_t2[s]         <- t2
    sam_xi[s]         <- xi
    sam_Theta[s,,]    <- Theta
    sam_s2[s]         <- s2
    sam_g[s,]         <- g
    sam_pg[s,]        <- pg
    sam_D[s,]         <- D
    
    ###########################################################################
    # Progress Bar Update
    ###########################################################################
    setTxtProgressBar(pb    = pb,
                      value = s / S)
    
  }
  
  # Returns
  return(list(B        = sam_B,
              l2       = sam_l2,
              l2_t2_s2 = sam_l2_t2_s2,
              l2_t2    = sam_l2_t2,
              e2_s2    = sam_e2_s2,
              v        = sam_v,
              t2       = sam_t2,
              xi       = sam_xi,
              Theta    = sam_Theta,
              s2       = sam_s2,
              g        = sam_g,
              pg       = sam_pg,
              D        = sam_D))
}