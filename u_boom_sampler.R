### Sampler

u_boom_sampler <- function(A,
                           G,
                           y,
                           S){
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
  m2       <- matrix(data = 1, nrow = P, ncol = P)
  diag(m2) <- 0
  w        <- matrix(data = 1, nrow = P, ncol = P)
  diag(w)  <- 0
  r2       <- 1
  nu       <- 1
  s2       <- 1
  
  # Sample Variables
  sam_B        <- array(data = NA, dim = c(S, V, P))
  sam_l2       <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2_s2 <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2    <- array(data = NA, dim = c(S, V, P))
  sam_v        <- array(data = NA, dim = c(S, V, P))
  sam_t2       <- numeric(length = S)
  sam_xi       <- numeric(length = S)
  sam_Theta    <- array(data = NA, dim = c(S, P, P))
  sam_m2       <- array(data = NA, dim = c(S, P, P))
  sam_m2_r2_s2 <- array(data = NA, dim = c(S, P, P))
  sam_m2_r2    <- array(data = NA, dim = c(S, P, P))
  sam_w        <- array(data = NA, dim = c(S, P, P))
  sam_r2       <- numeric(length = S)
  sam_nu       <- numeric(length = S)
  sam_s2       <- numeric(length = S)
  sam_g        <- matrix(data = NA, nrow = S, ncol = P)
  sam_pg       <- matrix(data = NA, nrow = S, ncol = P)
  
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
    out <- u_boom_iteration(y  = y,
                            X  = X,
                            l2 = l2,
                            v  = v,
                            t2 = t2,
                            xi = xi,
                            m2 = m2,
                            w  = w,
                            r2 = r2,
                            nu = nu,
                            s2 = s2,
                            g  = rep(1, P),
                            e2 = 1,
                            o2 = 1)
    # Assigns and Distributes the values
    # B Structure
    B     <- out$B
    l2    <- out$l2
    v     <- out$v
    t2    <- out$t2
    xi    <- out$xi
    # Theta Structure
    Theta <- out$Theta
    m2    <- out$m2
    w     <- out$w
    r2    <- out$r2
    nu    <- out$nu
    # Error
    s2 <- out$s2
    # Sparsity Structure
    g  <- out$g
    pg <- out$pg
    
    # Saves Variables
    sam_B[s,,]        <- B
    sam_l2[s,,]       <- l2
    sam_l2_t2_s2[s,,] <- l2 * t2 * s2
    sam_l2_t2[s,,]    <- l2 * t2
    sam_v[s,,]        <- v
    sam_t2[s]         <- t2
    sam_xi[s]         <- xi
    sam_Theta[s,,]    <- Theta
    sam_m2[s,,]       <- m2
    sam_m2_r2_s2[s,,] <- m2 * r2 * s2
    sam_m2_r2[s,,]    <- m2 * r2
    sam_w[s,,]        <- w
    sam_r2[s]         <- r2
    sam_nu[s]         <- nu
    sam_s2[s]         <- s2
    sam_g[s,]         <- g
    sam_pg[s,]        <- pg
    
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
              v        = sam_v,
              t2       = sam_t2,
              xi       = sam_xi,
              Theta    = sam_Theta,
              m2       = sam_m2,
              m2_r2_s2 = sam_m2_r2_s2,
              m2_r2    = sam_m2_r2,
              w        = sam_w,
              r2       = sam_r2,
              nu       = sam_nu,
              s2       = sam_s2,
              g        = sam_g,
              pg       = sam_pg))
}