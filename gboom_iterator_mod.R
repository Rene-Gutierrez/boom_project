gboom_iterator_mod <- function(y,
                               GX,
                               AX,
                               B,
                               C,
                               Theta,
                               s2 = 1,
                               l2B,
                               r2B,
                               t2B,
                               vB,
                               wB,
                               xiB,
                               l2T,
                               t2T,
                               vT,
                               xiT,
                               g,
                               ar1,
                               ar2,
                               r  = FALSE,
                               as,
                               tM,
                               bM,
                               full,
                               samr){
  # Problem dimensions
  n <- length(y)  # Number of Observations
  V <- dim(B)[1]  # Voxel Size
  P <- dim(B)[2]  # Number of ROI's
  
  # Samples r
  if(samr){
    r <- rbeta(n      = 1,
               shape1 = sum(g) + ar1,
               shape2 = P - sum(g) + ar2)
  } else {
    r <- 1 / 2
  }
  
  # Samples g (along with B and Theta)
  # Auxiliary Variables
  if(full){
    LB  <- t(t(l2B) * r2B) * t2B
    LT  <- l2T * t2T
    per <- sample(1:P)
    for(p in per){
      out <- group_iterator_mod(AX    = AX,
                                GX    = GX,
                                y     = y,
                                Theta = Theta,
                                LT    = LT,
                                B     = B,
                                C     = C,
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
  } else {
    g <- rep(1, P)
  }
  
  # Auxiliary Variables
  ipb <- bM[, g == 1][C[, g == 1]]
  ipt <- tM[g == 1, g == 1]
  ipt <- ipt[lower.tri(ipt)]
  iqb <- bM[, g != 1][C[, g != 1]]
  iqt <- seq(1, P * (P - 1) / 2)[-ipt]
  
  # Samples B and Theta jointly given g
  X <- cbind(GX[, ipb], AX[, ipt])
  temp <- t(r2B * t(l2B)) * t2B
  L <- c(temp[ipb], t2T * l2T[lower.tri(l2T)][ipt])
  # Jointly Fast Sampling of B and Theta
  b <- tryCatch(
    {
      fast_sampler(Phi = X / sqrt(s2),
                   D   = s2 * L,
                   a   = y / sqrt(s2))
    },
    error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      return(NA)
    }
  )
  if(!is.na(b[1])){
    # Updates B and Theta
    B[ ,g == 1][C[ ,g == 1]]     <- b[1:length(ipb)]
    Theta[upper.tri(Theta)]      <- 0
    Theta[lower.tri(Theta)][ipt] <- b[(length(ipb) + 1): length(c(ipb, ipt))]
    Theta <- Theta + t(Theta)
  }
  
  # Horseshoe Structure for B
  # Samples l2B
  for(p in 1:P){
    if(g[p] == 1){
      l2B[C[, p], p] <- 1 / rgamma(n     = sum(C[, p]),
                                   shape = 1,
                                   rate  = 1 / vB[C[, p], p] + B[C[, p], p]^2 / (2 * s2 * r2B[p] * t2B))
    } else {
      l2B[C[, p], p] <- 1 / rgamma(n     = sum(C[, p]),
                             shape = 1 / 2,
                             rate  = 1 / vB[C[, p], p])
    }
  }
  
  # Samples r2B
  for(i in 1:P){
    if(g[i] > 0){
      r2B[i] <- 1 / rgamma(n     = 1,
                           shape = (sum(C[, i]) + 1) / 2,
                           rate  = 1 / wB[i] +
                             sum(B[C[, i], i]^2 / (2 * s2 * l2B[C[, i], i] * t2B)))
    } else {
      r2B[i] <- 1 / rgamma(n     = 1,
                           shape = 1 / 2,
                           rate  = 1 / wB[i])
    }
    r2B[i] <- 1
  }
  
  # Samples t2B
  if(sum(g) > 0){
    t2B <- 1 / rgamma(n     = 1,
                      shape = (sum(C[, g == 1]) + 1) / 2,
                      rate  = 1 / xiB +
                        sum(B[C[, g == 1]]^2 / (2 * s2 * t(t(l2B) * r2B)[C[, g == 1]])))
  } else {
    t2B <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiB)
  }
  
  # Samples v2B
  vB[C]  <- 1 / rgamma(n     = sum(C),
                       shape = 1,
                       rate  = 1 + 1 / l2B[C])
  # Samples w2B
  wB  <- 1 / rgamma(n     = P,
                    shape = 1,
                    rate  = 1 + 1 / r2B)
  
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
                        sum(Theta[upper.tri(Theta)]^2 / (2 * s2 * l2T[upper.tri(l2T)])))
  } else {
    t2T <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiT)
  }
  
  # Samples vT
  out <- 1 / rgamma(n     = P * (P - 1) / 2,
                    shape = 1,
                    rate  = 1 + 1 / l2T[lower.tri(l2T)])
  vT[upper.tri(vT)] <- 0
  vT[lower.tri(vT)] <- out
  vT                <- vT + t(vT)
  
  # Samples xiT
  xiT <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2T)
  
  # Samples s2
  # Auxiliary Variables
  Tv  <- Theta[lower.tri(Theta)]
  Bv  <- B[C]
  Q   <- sum(g)
  lBv <- t(r2B * t(l2B))[C] * t2B
  lTv <- t2T * l2T[lower.tri(l2T)]
  # Samples
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + sum(C) + Q * (Q - 1) / 2 + as) / 2,
                   rate  = crossprod(x = y - GX[,C] %*% Bv - AX %*% Tv) / 2 +
                     t(Bv / lBv) %*% Bv / 2 + t(Tv / lTv) %*% Tv / 2 + as / 2)
  return(list(Theta = Theta,
              B     = B,
              g     = g,
              r     = r,
              l2B   = l2B,
              r2B   = r2B,
              t2B   = t2B,
              vB    = vB,
              wB    = wB,
              xiB   = xiB,
              l2T   = l2T,
              t2T   = t2T,
              vT    = vT,
              xiT   = xiT,
              s2    = s2))
  
  # Debbuger return
  return(list(B     = B,
              l2B   = l2B,
              vB    = vB,
              r2B   = r2B,
              wB    = wB,
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
              iqt   = iqt,
              GX    = GX,
              AX    = AX,
              r     = r))
}