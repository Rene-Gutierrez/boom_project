# Data Generator
data_generator2 <- function(P  = 10,
                            V  = 100,
                            pB = 1/2,
                            pT = 1/2,
                            cB = 1,
                            cT = 1,
                            s2 = 1,
                            n  = 100){
  # Samples the Non-zero indicators
  gT <- rbinom(n = P, size = 1, prob = pT)
  if(sum(gT) == 0){
    gT[1] <- 1
    gT[2] <- 1
  } 
  if(sum(gT) == 1){
    if(gT[1] == 0){
      gT[1] <- 1
    } else {
      gT[2] <- 1
    }
  }
  # if(sum(gT) == 0){
  #   gT[1] <- 1
  # }
  gB <- matrix(data = NA, nrow = V, ncol = P)
  for(p in 1:P){
    gB[, p] <- gT[p] * rbinom(n = V, size = 1, prob = pB)
    if(sum(gB[, p]) == 0){
      gB[1, p] <- 1 * gT[p]
    }
  }
  
  
  # Samples Theta
  Theta <- (gT %*% t(gT)) * matrix(data = rnorm(n    = P * P,
                                                mean = cT,
                                                sd   = 1),
                                   nrow = P,
                                   ncol = P)
  Theta[upper.tri(Theta)] <- 0
  diag(Theta)             <- 0
  Theta                   <- Theta + t(Theta)
  
  # Samples B
  B <- gB * matrix(data = rnorm(n    = V * P,
                                mean = cB,
                                sd   = 1),
                   nrow = V,
                   ncol = P)
  
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
  G <- array(data = rnorm(n    = n * V * P,
                          mean = 0,
                          sd   = 1),
             dim  = c(n, V, P))
  
  # Response Variable
  y <- numeric(length = n)
  for(i in 1:n){
    y[i] <- sum(G[i,,] * B) + sum(A[i,,] * Theta) / 2 +
      rnorm(n = 1, mean = 0, sd = sqrt(s2))
  }
  
  # Number of Coefficients in the Model
  nC <- sum(Theta != 0) + sum(B != 0) 
  
  # Returns the Data Generated
  return(list(gT    = gT,
              gB    = gB,
              Theta = Theta,
              B     = B,
              A     = A,
              G     = G,
              y     = y,
              s2    = s2,
              nC    = nC))
}