sub_structure_iterator <- function(GX,
                                   AX,
                                   y,
                                   B,
                                   Theta,
                                   LB,
                                   s2,
                                   gT,
                                   gB,
                                   iB,
                                   rB,
                                   p){
  # Problem dimensions
  n  <- length(y)  # Number of Observations
  V  <- dim(B)     # Voxel Dimensions
  P  <- dim(B)[1]  # Number of ROI's
  V  <- V[-1]      # Voxel Dimension
  nV <- length(V)  # Number of Voxel Dimensions
  
  # Updates gB
  # Auxiliary Variables
  iT       <- gT %*% t(gT)
  diag(iT) <- 0
  auxV     <- rep(0, P)
  auxV[p]  <- 1
  auxA     <- array(data = 1, dim = dim(B))
  subB     <- (auxA * auxV) == 1
  # Selects Active Variables
  unsT <- iT[lower.tri(iT)]
  unsB <- iB * (auxA * (1 - auxV)) == 1
  
  # Rewrites the the Regression Equation
  if(sum(unsB) == 0){
    if(sum(unsT) == 0){
      S <- y
    } else if(sum(unsT) == 1){
      S <- y - AX[, unsT] * Theta[lower.tri(Theta)][unsT]
    } else {
      S <- y - AX[, unsT] %*% Theta[lower.tri(Theta)][unsT]
    }
  } else if(sum(unsB) == 1){
    if(sum(unsT) == 0){
      S <- y - GX[, c(unsB)] * B[unsB]
    } else if(sum(unsT) == 1){
      S <- y - GX[, c(unsB)] * B[unsB] - AX[, unsT] * Theta[lower.tri(Theta)][unsT]
    } else {
      S <- y - GX[, c(unsB)] * B[unsB] - AX[, unsT] %*% Theta[lower.tri(Theta)][unsT]
    }
  } else {
    if(sum(unsT) == 0){
      S <- y - GX[, c(unsB)] %*% c(B[unsB])
    } else if(sum(unsT) == 1){
      S <- y - GX[, c(unsB)] %*% c(B[unsB]) - AX[, unsT] * Theta[lower.tri(Theta)][unsT]
    } else {
      S <- y - GX[, c(unsB)] %*% c(B[unsB]) - AX[, unsT] %*% Theta[lower.tri(Theta)][unsT]
    }
  }
  
  Bpro <- list()
  Bod  <- list()
  Bbh  <- list()
  Bsel <- list()
  Buns <- list()
  for(i in 1:nV){
    nVV       <- length(gB[[p]][[i]])
    Bpro[[i]] <- numeric(length = nVV)
    Bod[[i]]  <- numeric(length = nVV)
    Bbh[[i]]  <- list()
    Bsel[[i]] <- list()
    Buns[[i]] <- list()
    for(j in 1:nVV){
      temp    <- gB[[p]][[i]]
      temp[j] <- 1
      if(i == 1){
        igB <- temp
      } else {
        igB <- gB[[p]][[1]]
      }
      if( nV > 1){
        for(k in 2:nV){
          if(k == i){
            igB <- igB %o% temp
          } else {
            igB <- igB %o% gB[[p]][[k]]
          }
        }
      }
      if(i == 1){
        sel    <- rep(0, nVV)
        sel[j] <- 1
      } else {
        sel <- rep(1, length(gB[[p]][[1]]))
      }
      if( nV > 1){
        for(k in 2:nV){
          if(k == i){
            tem    <- rep(0, nVV)
            tem[j] <- 1
            sel <- sel %o% tem
          } else {
            sel <- sel %o% rep(1, length(gB[[p]][[k]]))
          }
        }
      }
      # Selected and Unselected
      selB       <- (auxA *  auxV) == 1
      selB[selB] <- (igB * sel) == 1
      unsB       <- (auxA *  auxV) == 1
      unsB[unsB] <- (igB * (1 - sel)) == 1
      
      Bsel[[i]][[j]] <- selB
      Buns[[i]][[j]] <- unsB
      # Rewrites the the Regression Equation
      if(sum(unsB) == 0){
        R <- S
      } else if(sum(unsB) == 1){
        R <- S - GX[, c(unsB)] * B[unsB]
      } else {
        R <- S - GX[, c(unsB)] %*% c(B[unsB])
      }
      X <- GX[, c(selB)]
      
      # Number of Estimates
      q <- sum(selB)
      
      if(q == 0){
        bh <- c()
        Bbh[[i]][[j]] <- bh
        # Odds
        od          <- rB / (1 - rB)
        Bod[[i]][j] <- od
        # Probability
        pro          <- 1 / (1 + od)
        Bpro[[i]][j] <- pro
        # Updates gB
        gB[[p]][[i]][j] <- rbinom(n = 1, size = 1, prob = pro)
      } else if(q == 1){
        # Obtains the point estimate
        LL <- 1 / LB[selB]
        Q  <- c(t(X) %*% X + LL)
        U  <- sqrt(Q)
        bh <- c(t(X) %*% R / Q)
        Bbh[[i]][[j]] <- bh
        # Odds for gT == 1
        od <- - log(Q) / 2
        od <- od + t(bh) * Q * bh / s2 / 2
        od <- od - log(LL) / 2
        od <- od + log(rB) - log(1 - rB)
        od <- exp(od)
        Bod[[i]][j] <- od
        # Probability for gT == 1
        if(is.infinite(od)){
          pro <- 1
        } else {
          pro <- od / (1 + od)
        }
        Bpro[[i]][j] <- pro
        # Updates gB
        gB[[p]][[i]][j] <- rbinom(n = 1, size = 1, prob = pro)
        # Updates B and Theta
        if(gB[[p]][[i]][j] == 1){
          b            <- backsolve(U, rnorm(n = q)) + bh
          B[selB == 1] <- b
        } else {
          b            <- rep(0, q)
          B[selB == 1] <- b
        }
      } else {
        # Obtains the point estimate
        LL <- diag(1 / LB[selB])
        Q  <- t(X) %*% X + LL
        U  <- chol(Q)
        bh <- backsolve(U, forwardsolve(t(U), t(X) %*% R))
        Bbh[[i]][[j]] <- bh
        # Odds for gT == 1
        od <- - log(det(Q)) / 2
        od <- od + t(bh) %*% Q %*% bh / s2 / 2
        od <- od - log(prod(diag(LL))) / 2
        od <- od + log(rB) - log(1 - rB)
        od <- exp(od)
        Bod[[i]][j] <- od
        # Probability for gT == 1
        if(is.infinite(od)){
          pro <- 1
        } else {
          pro <- od / (1 + od)
        }
        Bpro[[i]][j] <- pro
        # Updates gB
        gB[[p]][[i]][j] <- rbinom(n = 1, size = 1, prob = pro)
        # Updates B and Theta
        if(gB[[p]][[i]][j] == 1){
          b            <- backsolve(U, rnorm(n = q)) + bh
          B[selB == 1] <- b
        } else {
          b            <- rep(0, q)
          B[selB == 1] <- b
        }
      }
    }
  }
  
  # Returns
  return(list(B    = B,
              Bod  = Bod,
              Bpro = Bpro,
              Bbh  = Bbh,
              Bsel = Bsel,
              Buns = Buns,
              gB   = gB))
}