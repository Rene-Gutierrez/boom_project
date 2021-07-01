pre_data_gen <- function(B,
                         Theta,
                         s2,
                         n = 100){
  # Dimensions
  V <- dim(B)[1]
  P <- dim(B)[2]
  
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
  
  # Returns the Data Generated
  return(list(A     = A,
              G     = G,
              y     = y))
}