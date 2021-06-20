gboom_iterator <- function(y,
                           GX,
                           AX,
                           B,
                           Theta,
                           s2,
                           l2B,
                           t2B,
                           vB,
                           xiB,
                           l2T,
                           t2T,
                           vT,
                           xiT,
                           g,
                           r,
                           tM,
                           bM){
  # Problem dimensions
  n <- length(y)  # Number of Observations
  V <- dim(B)[1]  # Voxel Size
  P <- dim(B)[2]  # Number of ROI's
  Q <- sum(g)     # Number of Active Regions
  
  # Samples g (along with B and Theta)
  # Auxiliary Variables
  LB <- l2B * t2B
  LT <- l2T * t2T
  per <- sample(1:P)
  for(p in per){
    out <- group_iterator(AX    = AX,
                          GX    = GX,
                          y     = y,
                          Theta = Theta,
                          LT    = LT,
                          B     = B,
                          LB    = LB,
                          s2    = s2,
                          g     = g,
                          r     = r,
                          p     = p,
                          tM    = tM,
                          bM    = bM)
    Theta <- out$Theta
    B     <- out$B
    g     <- out$g
  }
  
  # Auxiliary Variables
  ipb <- c(bM[, g == 1])
  ipt <- tM[g == 1, g == 1]
  ipt <- ipt[lower.tri(ipt)]
  iqb <- c(bM[, g != 1])
  iqt <- seq(1, P * (P - 1) / 2)[-ipt]
  
  # Samples B and Theta jointly given g
  X <- cbind(GX[, ipb], AX[, ipt])
  L <- c(t2B * c(l2B)[ipb], t2T * l2T[lower.tri(l2T)][ipt])
  # Jointly Fast Sampling of B and Theta
  b <- fast_sampler(Phi = X / sqrt(s2),
                    D   = s2 * L,
                    a   = y / sqrt(s2))
  # Updates B and Theta
  B[ ,g == 1]                  <- b[1:length(ipb)]
  Theta[upper.tri(Theta)]      <- 0
  Theta[lower.tri(Theta)][ipt] <- b[(length(ipb) + 1): length(c(ipb, ipt))]
  Theta <- Theta + t(Theta)
  
  # Horseshoe Structure for B
  # Samples l2B
  for(i in 1:P){
    if(g[i] == 1){
      l2B[, i] <- 1 / rgamma(n     = V,
                             shape = 1,
                             rate  = 1 / vB[, i] + B[, i]^2 / (2 * s2 * t2B))
    } else {
      l2B[, i] <- 1 / rgamma(n     = V,
                             shape = 1 / 2,
                             rate  = 1 / vB[, i])
    }
  }
  
  # Samples t2B
  if(sum(g) > 0){
    t2B <- 1 / rgamma(n     = 1,
                      shape = (sum(g) * V + 1) / 2,
                      rate  = 1 / xiB +
                        sum(B^2 / (2 * s2 * l2B)))
  } else {
    t2B <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiB)
  }
  
  # Samples v2B
  vB  <- 1 / rgamma(n     = P * V,
                    shape = 1,
                    rate  = 1 + 1 / c(l2B))
  vB <- matrix(data = vB, nrow = V, ncol = P)
  
  # Samples xiB
  xiB <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2B)
  
  # Horseshoe Structure for Theta
  # Samples l2T
  for(i in 2:P){
    for(j in 1:(i - 1)){
      if(g[i] * g[j] == 1){
        l2T[i, j] <- 1 / rgamma(n     = 1,
                                shape = 1,
                                rate  = 1 / vT[i, j] + Theta[i, j]^2 / (2 * s2 * t2T))
        l2T[j, i] <- l2T[i, j]
      } else {
        l2T[i, j] <- 1 / rgamma(n     = 1,
                                shape = 1 / 2,
                                rate  = 1 / vT[i, j])
        l2T[j, i] <- l2T[i, j]
      }
    }
  }
  
  # Samples t2T
  if(sum(g) > 0){
    t2T <- 1 / rgamma(n     = 1,
                      shape = (sum(g) * (sum(g) - 1) / 2 + 1) / 2,
                      rate  = 1 / xiT +
                        sum(Theta^2 / (2 * s2 * l2T)))
  } else {
    t2T <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiT)
  }
  
  # Samples vT
  out <- 1 / rgamma(n     = P * (P - 1) / 2,
                    shape = 1,
                    rate  = 1 + 1 / l2B[lower.tri(l2B)])
  vT[upper.tri(vT, diag = TRUE)] <- 0
  vT[lower.tri(vT)]              <- out
  vT                             <- vT + t(vT)
  
  # Samples xiT
  xiT <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2T)
  
  # Samples s2
  # Auxilary Variables
  Tv  <- Theta[lower.tri(Theta)]
  Bv  <- c(B)
  Q   <- sum(g)
  lBv <- t2B * c(l2B)
  lTv <- t2T * l2T[lower.tri(l2T)]
  # Samples
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + Q * V + Q * (Q - 1) / 2) / 2,
                   rate  = crossprod(x = y - GX %*% Bv - AX %*% Tv) / 2 +
                     t(Bv / lBv) %*% Bv / 2 + t(Tv / lTv) %*% Tv / 2)
  
  
  return(list(B     = B,
              l2B   = l2B,
              vB    = vB,
              t2B   = t2B,
              xiB   = xiB,
              Theta = Theta,
              l2T   = l2T,
              vT    = vT,
              t2T   = t2T,
              xiT   = xiT,
              s2    = s2,
              g     = g,
              ipb   = ipb,
              ipt   = ipt,
              iqb   = iqb,
              iqt   = iqt))
}