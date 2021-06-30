full_group_iterator <- function(AX,
                                GX,
                                y,
                                Theta,
                                LT,
                                B,
                                LB,
                                s2,
                                gT,
                                gB,
                                rT,
                                rB,
                                p){
  # Problem dimensions
  n  <- length(y)  # Number of Observations
  V  <- dim(B)     # Voxel Dimensions
  P  <- dim(B)[1]  # Number of ROI's
  V  <- V[-1]      # Voxel Dimension
  nV <- length(V)  # Number of Voxel Dimensions
  
  # Auxiliary Variables
  gTemp    <- gT
  gTemp[p] <- 1
  iT       <- gTemp %*% t(gTemp)
  diag(iT) <- 0
  iB       <- array(data = NA, dim = c(0, V))
  for(i in 1:P){
    igB <- gB[[i]][[1]]
    if( nV > 1){
      for(j in 2:nV){
        igB <- igB %o% gB[[i]][[j]] 
      }
    }
    igB <- igB * gTemp[i]
    iB  <- abind::abind(iB, igB, along = 1)
  }
  
  # gT Sampling
  auxp    <- rep(1, P)
  auxp[p] <- 0
  # Selected from Theta
  selT <- t(iT * (1 - auxp %*% t(auxp)))[lower.tri(iT)] == 1
  unsT <- t(iT * (auxp %*% t(auxp)))[lower.tri(iT)] == 1
  # Selections from B
  selB <- (iB * (1 - auxp)) == 1
  unsB <- (iB * auxp) == 1
  
  # Rewrites the the Regression Equation
  if(sum(unsB) == 0){
    if(sum(unsT) == 0){
      R <- y
    } else if(sum(unsT) == 1){
      R <- y - AX[, unsT] * Theta[lower.tri(Theta)][unsT]
    } else {
      R <- y - AX[, unsT] %*% Theta[lower.tri(Theta)][unsT]
    }
  } else if(sum(unsB) == 1){
    if(sum(unsT) == 0){
      R <- y - GX[, c(unsB)] * B[unsB]
    } else if(sum(unsT) == 1){
      R <- y - GX[, c(unsB)] * B[unsB] - AX[, unsT] * Theta[lower.tri(Theta)][unsT]
    } else {
      R <- y - GX[, c(unsB)] * B[unsB] - AX[, unsT] %*% Theta[lower.tri(Theta)][unsT]
    }
  } else {
    if(sum(unsT) == 0){
      R <- y - GX[, c(unsB)] %*% c(B[unsB])
    } else if(sum(unsT) == 1){
      R <- y - GX[, c(unsB)] %*% c(B[unsB]) - AX[, unsT] * Theta[lower.tri(Theta)][unsT]
    } else {
      R <- y - GX[, c(unsB)] %*% c(B[unsB]) - AX[, unsT] %*% Theta[lower.tri(Theta)][unsT]
    }
  }
  X <- cbind(GX[, c(selB)], AX[, selT])
  
  # Number of Estimates
  q <- sum(selB) + sum(selT)
  
  if(q == 0){
    bh <- c()
    # Odds
    od  <- rT / (1 - rT)
    # Probability
    pro <- 1 / (1 + od)
    # Updates gT
    gT[p] <- rbinom(n = 1, size = 1, prob = pro)
  } else if(q == 1){
    # Obtains the point estimate
    LL <- 1 / c(LB[selB], LT[lower.tri(LT)][selT])
    Q  <- c(t(X) %*% X + LL)
    U  <- sqrt(Q)
    bh <- c(t(X) %*% R / Q)
    # Odds for gT == 1
    od <- - log(Q) / 2
    od <- od + t(bh) * Q * bh / s2 / 2
    od <- od - log(LL) / 2
    od <- od + log(rT) - log(1 - rT)
    od <- exp(od)
    # Probability for gT == 1
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    # Updates gT
    gT[p] <- rbinom(n = 1, size = 1, prob = pro)
    # Updates B and Theta
    if(gT[p] == 1){
      b <- backsolve(U, rnorm(n = q)) + bh
      if(sum(selB) > 0){
        B[selB] <- b[1:sum(selB)]
      }
      if(sum(selT) > 0){
        Theta[lower.tri(Theta)][selT]        <- b[(sum(selB) + 1):q]
        Theta[upper.tri(Theta, diag = TRUE)] <- 0
        Theta <- Theta + t(Theta)
      }
    } else {
      b <- rep(0, q)
      if(sum(selB) > 0){
        B[selB] <- b[1:sum(selB)]
      }
      if(sum(selT) > 0){
        Theta[lower.tri(Theta)][selT]        <- b[(sum(selB) + 1):q]
        Theta[upper.tri(Theta, diag = TRUE)] <- 0
        Theta <- Theta + t(Theta)
      }
    }
  } else {
    # Obtains the point estimate
    LL <- diag(1 / c(LB[selB], LT[lower.tri(LT)][selT]))
    Q  <- t(X) %*% X + LL
    U  <- chol(Q)
    bh <- backsolve(U, forwardsolve(t(U), t(X) %*% R))
    # Odds for gT == 1
    od <- - log(det(Q)) / 2
    od <- od + t(bh) %*% Q %*% bh / s2 / 2
    od <- od - log(prod(diag(LL))) / 2
    od <- od + log(rT) - log(1 - rT)
    od <- exp(od)
    # Probability for gT == 1
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    # Updates g
    gT[p] <- rbinom(n = 1, size = 1, prob = pro)
    # Updates B and Theta
    if(gT[p] == 1){
      b <- backsolve(U, rnorm(n = q)) + bh
      if(sum(selB) > 0){
        B[selB] <- b[1:sum(selB)]
      }
      if(sum(selT) > 0){
        Theta[lower.tri(Theta)][selT]        <- b[(sum(selB) + 1):q]
        Theta[upper.tri(Theta, diag = TRUE)] <- 0
        Theta <- Theta + t(Theta)
      }
    } else {
      b <- rep(0, q)
      if(sum(selB) > 0){
        B[selB] <- b[1:sum(selB)]
      }
      if(sum(selT) > 0){
        Theta[lower.tri(Theta)][selT]        <- b[(sum(selB) + 1):q]
        Theta[upper.tri(Theta, diag = TRUE)] <- 0
        Theta <- Theta + t(Theta)
      }
    }
  }
  
  if(gT[p] == 1){
    Bsub <- sub_structure_iterator(GX    = GX,
                                   AX    = AX,
                                   y     = y,
                                   B     = B,
                                   Theta = Theta,
                                   LB    = LB,
                                   s2    = s2,
                                   gT    = gT,
                                   gB    = gB,
                                   iB    = iB,
                                   rB    = rB,
                                   p     = p)
    TB  <- B
    B   <- Bsub$B
    TgB <- gB
    gB  <- Bsub$gB
  } else {
    Bsub <- list()
    TB   <- B
    TgB  <- gB
    for(i in 1:nV){
      nVV <- length(gB[[i]])
      if(nVV > 1){
        for(j in 1:nVV){
          gB[[p]][[i]][j] <- rbinom(n = 1, size = 1, prob = rB)
        }
      }
    }
  }
  
  return(list(n     = n,
              V     = V,
              P     = P,
              iT    = iT,
              iB    = iB,
              selT  = selT,
              unsT  = unsT,
              selB  = selB,
              unsB  = unsB,
              R     = R,
              X     = X,
              pro   = pro,
              od    = od,
              bh    = bh,
              Theta = Theta,
              B     = B,
              TB    = TB,
              TgB   = TgB,
              gT    = gT,
              gB    = gB,
              sub   = Bsub))
}