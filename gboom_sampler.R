gboom_sampler <- function(y,
                          G,
                          A,
                          g,
                          s2,
                          r,
                          S){
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
  sam_B     <- array(data = NA, dim = c(S, V, P))
  sam_Theta <- array(data = NA, dim = c(S, P, P))
  sam_s2    <- numeric(length = S)
  sam_l2B   <- array(data = NA, dim = c(S, V, P))
  sam_t2B   <- numeric(length = S)
  sam_vB    <- array(data = NA, dim = c(S, V, P))
  sam_xiB   <- numeric(length = S)
  sam_l2T   <- array(data = NA, dim = c(S, P, P))
  sam_t2T   <- numeric(length = S)
  sam_vT    <- array(data = NA, dim = c(S, P, P))
  sam_xiT   <- numeric(length = S)
  sam_g     <- matrix(data = NA, nrow = S, ncol = P)
  
  # Variable Initialization
  B     <- matrix(data = 0, nrow = V, ncol = P)
  Theta <- matrix(data = 0, nrow = P, ncol = P)
  l2B   <- matrix(data = 1, nrow = V, ncol = P)
  t2B   <- 1
  vB    <- matrix(data = 1, nrow = V, ncol = P)
  xiB   <- 1
  l2T   <- matrix(data = 1, nrow = P, ncol = P)
  t2T   <- 1
  vT    <- matrix(data = 1, nrow = P, ncol = P)
  xiT   <- 1
  
  # Sampling
  # Progress Bar
  pb <- txtProgressBar(min     = 0,
                       max     = 1,
                       initial = 0,
                       style   = 3,
                       width   = 72)
  for(s in 1:S){
    # Samples
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
                          bM     = bM)
    
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
    sam_B[s,,]     <- B
    sam_Theta[s,,] <- Theta
    sam_s2[s]      <- s2
    sam_l2B[s,,]   <- l2B
    sam_t2B[s]     <- t2B
    sam_vB[s,,]    <- vB
    sam_xiB[s]     <- xiT
    sam_l2T[s,,]   <- l2T
    sam_t2T[s]     <- t2T
    sam_vT[s,,]    <- vT
    sam_xiT[s]     <- xiT
    sam_g[s,]      <- g
    
    # Progress Bar Update
    setTxtProgressBar(pb    = pb,
                      value = s / S)
  }
  
  # Returns the Samples
  return(list(B     = sam_B,
              Theta = sam_Theta,
              s2    = sam_s2,
              g     = sam_g,
              l2B   = sam_l2B,
              t2B   = sam_t2B,
              vB    = sam_vB,
              xiB   = sam_xiB,
              l2T   = sam_l2T,
              t2T   = sam_t2T,
              vT    = sam_vT,
              xiT   = sam_xiT))
}