### Sampler

horseshoe_sampler <- function(A,
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
  GX                <- G
  dim(GX)           <- c(n, V * P)
  AX                <- A
  dim(AX)           <- c(n, P * P)
  AX                <- AX[, lower.tri(A[1,,])]
  X                 <- cbind(GX, AX)
  y                 <- y
  
  # Variable Initialization
  xi <- 1
  v  <- rep(1, P * V + P * (P - 1) / 2)
  t2 <- 1
  l2 <- rep(1, P * V + P * (P - 1) / 2)
  s2 <- 1
  
  # Sample Variables
  sam_B         <- array(data = NA, dim = c(S, V, P))
  sam_l2B       <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2_s2B <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2B    <- array(data = NA, dim = c(S, V, P))
  sam_vB        <- array(data = NA, dim = c(S, V, P))
  sam_Theta     <- array(data = NA, dim = c(S, P, P))
  sam_l2T       <- array(data = NA, dim = c(S, P, P))
  sam_l2_t2_s2T <- array(data = NA, dim = c(S, P, P))
  sam_l2_t2T    <- array(data = NA, dim = c(S, P, P))
  sam_vT        <- array(data = NA, dim = c(S, P, P))
  sam_t2        <- numeric(length = S)
  sam_xi        <- numeric(length = S)
  sam_s2        <- numeric(length = S)
  
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
    out <- fast_horseshoe(y  = y,
                          X  = X,
                          xi = xi,
                          v  = v,
                          l2 = l2,
                          t2 = t2,
                          s2 = s2)
    # Retrieves the Values
    b  <- out$b
    l2 <- out$l2
    t2 <- out$t2
    s2 <- out$s2
    v  <- out$v
    xi <- out$xi
    
    # Assigns the values
    B     <- matrix(data = out$b[1:(V * P)], V, P)
    Theta <- matrix(data = 0, nrow = P, ncol = P)
    temp  <- out$b[(V * P + 1):(V * P + (P - 1) * P / 2)]
    Theta[lower.tri(Theta)] <- temp
    Theta <- Theta + t(Theta)
    s2    <- out$s2
    l2B   <- matrix(data = out$l2[1:(V * P)], V, P)
    l2T   <- matrix(data = 0, nrow = P, ncol = P)
    temp  <- out$l2[(V * P + 1):(V * P + (P - 1) * P / 2)]
    l2T[lower.tri(l2T)] <- temp
    l2T   <- l2T + t(l2T)
    vB    <- matrix(data = out$v[1:(V * P)], V, P)
    vT    <- matrix(data = 0, nrow = P, ncol = P)
    temp  <- out$v[(V * P + 1):(V * P + (P - 1) * P / 2)]
    vT[lower.tri(vT)] <- temp
    vT    <- vT + t(vT)
    xi    <- out$xi
    t2    <- out$t2
    
    # Saves Variables
    sam_B[s,,]         <- B
    sam_Theta[s,,]     <- Theta
    sam_l2B[s,,]       <- l2B
    sam_vB[s,,]        <- vB
    sam_l2_t2B[s,,]    <- l2B * t2
    sam_l2_t2_s2B[s,,] <- l2B * t2 * s2
    sam_l2T[s,,]       <- l2T
    sam_vT[s,,]        <- vT
    sam_l2_t2T[s,,]    <- l2T * t2
    sam_l2_t2_s2T[s,,] <- l2T * t2 * s2
    sam_t2[s]          <- t2
    sam_xi[s]          <- xi
    sam_s2[s]          <- s2
    
    ###########################################################################
    # Progress Bar Update
    ###########################################################################
    setTxtProgressBar(pb    = pb,
                      value = s / S)
    
  }
  
  # Returns
  return(list(B         = sam_B,
              Theta     = sam_Theta,
              s2        = sam_s2,
              l2B       = sam_l2B,
              l2T       = sam_l2T,
              t2        = sam_t2,
              vB        = sam_vB,
              vT        = sam_vT,
              xi        = sam_xi,
              l2_t2B    = sam_l2_t2B,
              l2_t2_s2B = sam_l2_t2_s2B,
              l2_t2T    = sam_l2_t2T,
              l2_t2_s2T = sam_l2_t2_s2T))
}