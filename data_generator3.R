# Data Generator
data_generator3 <- function(P  = 10,
                            V  = c(4,5),
                            pB = 1/2,
                            pT = 1/2,
                            pT2 = 1/2,
                            cB = 1,
                            cT = 1,
                            s2 = 1,
                            n  = 100){
  # Activation Indicators
  iT <- matrix(data = 0, nrow = P, ncol = P) 
  # Samples the Non-zero indicators
  gT <- rbinom(n = P, size = 1, prob = pT)
  for(p in 2:P){
    for(q in 1:(p - 1)){
      if((gT[p] * gT[q]) == 1){
        iT[p, q] <- rbinom(n = 1, size = 1, prob = pT2)
      }
    }
  }
  
  # Samples Theta
  Theta <- iT * matrix(data = rnorm(n    = P * P,
                                    mean = cT,
                                    sd   = 1),
                       nrow = P,
                       ncol = P)
  # Makes the Matrices Symmetric
  Theta <- Theta + t(Theta)
  iT    <- iT + t(iT)
  
  # Recomputes gT based on real sparsity level
  for(p in 1:P){
    if(sum(iT[, p]) > 0){
      gT[p] <- 1
    } else {
      gT[p] <- 0
    }
  }
  
  iB <- array(data = NA, dim = c(0, V))
  gB <- list()
  for(p in 1:P){
    nV      <- length(V)
    gB[[p]] <- list()
    for(i in 1:nV){
      gB[[p]][[i]] <- rbinom(n = V[i], size = 1, prob = pB)
    }
    igB <- gB[[p]][[1]]
    if( nV > 1){
      for(i in 2:nV){
        igB <- igB %o% gB[[p]][[i]] 
      }
    }
    if(gT[p] == 0){
      igB <- igB * 0
    }
    iB <- abind::abind(iB, igB, along = 1)
  }
  
  # Samples B
  B <- iB * array(data = rnorm(n    = P * prod(V), 
                               mean = cB,
                               sd   = 1),
                  dim  = c(P, V))
  
  # Covariates
  # A
  A <- array(data = rnorm(n    = n * P * P,
                          mean = 0,
                          sd   = 1 / 2),
             dim  = c(n, P, P))
  for(i in 1:n){
    diag(A[i,,]) <- 0
    A[i,,]       <- A[i,,] + t(A[i,,])
  }
  # G
  G <- array(data = rnorm(n    = n * P * prod(V),
                          mean = 0,
                          sd   = 1),
             dim  = c(n, P, V))
  
  # Response Variable
  y <- apply(X = G, MARGIN = 1, FUN = function(x) sum(x * B))
  for(i in 1:n){
    y[i] <- y[i] + sum(A[i,,] * Theta) / 2
    y[i] <- y[i] + rnorm(n = 1, mean = 0, sd = sqrt(s2))
  }
  
  # Number of Coefficients in the Model
  nC <- sum(Theta != 0) + sum(B != 0) 
  
  # Returns the Data Generated
  return(list(gT    = gT,
              iT    = iT,
              gB    = gB,
              iB    = iB,
              Theta = Theta,
              B     = B,
              A     = A,
              G     = G,
              y     = y,
              s2    = s2,
              nC    = nC))
}