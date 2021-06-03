group_theta_sampler <- function(A, y, g, Theta, L, s, r, S){
  # Dimensions
  P <- dim(Theta)[1]
  n <- length(y)
  X <- A
  dim(X) <- c(n, P * P)
  X <- X[, lower.tri(Theta)]
  
  # Initialization
  for(p in 1:P){
    if(g[p] == 0){
      Theta[p, ] <- 0
      Theta[, p] <- 0
    }
  }
  # Samples
  sam_Theta <- array(data =NA, dim = c(S, P, P))
  sam_g     <- matrix(data = NA, nrow = S, ncol = P)
  sam_s     <- numeric(length = S)
  
  for(i in 1:S){
    per <- sample(1:P)
    for(p in per){
      print(p)
      out <- groupSymSpikeSlabUpd(A     = A,
                                  y     = y,
                                  Theta = Theta,
                                  L     = L,
                                  s     = s,
                                  g     = g,
                                  p     = p,
                                  r     = r)
      Theta <- out$Theta
      g     <- out$g
      print(out$pro)
      
      # Updates s
      s <- 1 / rgamma(n     = 1,
                      shape = (n + sum(g[p]) * (sum(g[p]) - 1) / 2) / 2,
                      rate  = crossprod(x = y - X %*% Theta[lower.tri(Theta)]) / 2 +
                        sum((Theta * Theta) / L))
      
    }
    sam_Theta[i,,] <- Theta
    sam_s[i]       <- s
    sam_g[i,]      <- out$g
  }
  return(list(Theta = sam_Theta,
              s     = sam_s,
              g     = sam_g))
}