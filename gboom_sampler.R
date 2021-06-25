gboom_sampler <- function(y,
                          G,
                          A,
                          g,
                          s2,
                          r,
                          S,
                          R = 0){
  # Problem Dimensions
  n <- dim(G)[1]  # Number of Observations
  V <- dim(G)[2]  # Voxel Size
  P <- dim(G)[3]  # Number of ROI's
  Q <- sum(g)     # Number of Active Regions
  
  # Auxilary Variables
  GX                <- G
  dim(GX)           <- c(nn, V * P)
  AX                <- A
  dim(AX)           <- c(n, P * P)
  AX                <- AX[, lower.tri(A[1,,])]
  tM                <- matrix(data = 0, nrow = P, ncol = P)
  tM[lower.tri(tM)] <- seq(1, (P - 1) * P / 2)
  tM                <- tM + t(tM)
  bM                <- matrix(data = 1:(P * V), nrow = V, ncol = P)
  
  # Samples
  sam_B         <- array(data = NA, dim = c(S, V, P))
  sam_Theta     <- array(data = NA, dim = c(S, P, P))
  sam_s2        <- numeric(length = S)
  sam_l2B       <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2_s2B <- array(data = NA, dim = c(S, V, P))
  sam_l2_t2B    <- array(data = NA, dim = c(S, V, P))
  sam_t2B       <- matrix(data = NA, nrow = S, ncol = P)
  sam_vB        <- array(data = NA, dim = c(S, V, P))
  sam_xiB       <- matrix(data = NA, nrow = S, ncol = P)
  sam_l2T       <- array(data = NA, dim = c(S, P, P))
  sam_l2_t2_s2T <- array(data = NA, dim = c(S, P, P))
  sam_l2_t2T    <- array(data = NA, dim = c(S, P, P))
  sam_t2T       <- numeric(length = S)
  sam_vT        <- array(data = NA, dim = c(S, P, P))
  sam_xiT       <- numeric(length = S)
  sam_g         <- matrix(data = NA, nrow = S, ncol = P)
  
  # Variable Initialization
  B     <- matrix(data = 0, nrow = V, ncol = P)
  Theta <- matrix(data = 0, nrow = P, ncol = P)
  l2B   <- matrix(data = 1, nrow = V, ncol = P)
  t2B   <- rep(1, P)
  vB    <- matrix(data = 1, nrow = V, ncol = P)
  xiB   <- rep(1, P)
  l2T   <- matrix(data = 1, nrow = P, ncol = P)
  t2T   <- 1
  vT    <- matrix(data = 1, nrow = P, ncol = P)
  xiT   <- 1
  iniS2 <- s2
  lg    <- g
  
  # Sampling
  # Progress Bar
  pb <- txtProgressBar(min     = 0,
                       max     = 1,
                       initial = 0,
                       style   = 3,
                       width   = 72)
  for(s in 1:S){
    # Samples
    if(s < R){
      full <- FALSE
    } else {
      full <- TRUE
    }
    out <- gboom_iterator(y      = y,
                          GX     = GX,
                          AX     = AX,
                          B      = B,
                          Theta  = Theta,
                          s2     = s2,
                          l2B    = l2B,
                          t2B    = t2B,
                          vB     = vB,
                          xiB    = xiB,
                          l2T    = l2T,
                          t2T    = t2T,
                          vT     = vT,
                          xiT    = xiT,
                          g      = g,
                          r      = r,
                          tM     = tM,
                          bM     = bM,
                          full   = full)
    
    # Updates the Values
    B     <- out$B
    Theta <- out$Theta
    s2    <- out$s2
    l2B   <- out$l2B
    t2B   <- out$t2B
    vB    <- out$vB
    xiB   <- out$xiB
    l2T   <- out$l2T
    t2T   <- out$t2T
    vT    <- out$vT
    xiT   <- out$xiT
    g     <- out$g
    
    # Saves the Samples
    sam_B[s,,]         <- B
    sam_Theta[s,,]     <- Theta
    sam_s2[s]          <- s2
    sam_l2B[s,,]       <- l2B
    sam_l2_t2_s2B[s,,] <- t(t(l2B) * t2B) * s2
    sam_l2_t2B[s,,]    <- t(t(l2B) * t2B)
    sam_t2B[s,]        <- t2B
    sam_vB[s,,]        <- vB
    sam_xiB[s]         <- xiT
    sam_l2T[s,,]       <- l2T
    sam_l2_t2_s2T[s,,] <- l2T * t2T * s2
    sam_l2_t2T[s,,]    <- l2T * t2T
    sam_t2T[s]         <- t2T
    sam_vT[s,,]        <- vT
    sam_xiT[s]         <- xiT
    sam_g[s,]          <- g
    
    # Re-Initialize
    if(sum(g) == 0){
      B     <- matrix(data = 0, nrow = V, ncol = P)
      Theta <- matrix(data = 0, nrow = P, ncol = P)
      l2B   <- matrix(data = 1, nrow = V, ncol = P)
      t2B   <- rep(1, P)
      vB    <- matrix(data = 1, nrow = V, ncol = P)
      xiB   <- rep(1, P)
      l2T   <- matrix(data = 1, nrow = P, ncol = P)
      t2T   <- 1
      vT    <- matrix(data = 1, nrow = P, ncol = P)
      xiT   <- 1
      s2    <- iniS2
      g     <- rbinom(n = P, size = 1, prob = 1 / 10)
    }
    
    # Progress Bar Update
    setTxtProgressBar(pb    = pb,
                      value = s / S)
  }
  
  # Returns the Samples
  return(list(B         = sam_B,
              Theta     = sam_Theta,
              s2        = sam_s2,
              g         = sam_g,
              l2B       = sam_l2B,
              l2_t2_s2B = sam_l2_t2_s2B,
              l2_t2B    = sam_l2_t2B,
              t2B       = sam_t2B,
              vB        = sam_vB,
              xiB       = sam_xiB,
              l2T       = sam_l2T,
              l2_t2_s2T = sam_l2_t2_s2T,
              l2_t2T    = sam_l2_t2T,
              t2T       = sam_t2T,
              vT        = sam_vT,
              xiT       = sam_xiT))
}