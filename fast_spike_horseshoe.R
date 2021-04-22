# Fast Horse Shoe Sampler (For Big p and small n)

fast_spike_horseshoe <- function(y,
                                 X,
                                 xi,
                                 v,
                                 l2,
                                 t2,
                                 s2,
                                 g,
                                 e2){
  #############################################################################
  # Dimensions
  #############################################################################
  n <- dim(X)[1]
  p <- dim(X)[2]
  P <- length(g)
  V <- p / P
  
  #############################################################################
  # Samples b (Fast)
  #############################################################################
  G <- c(t(matrix(data = rep(g, V), nrow = P, ncol = V)))
  t <- c(t(matrix(data = rep(t2, V), nrow = P, ncol = V)))
  E <- c(t(matrix(data = rep(e2, V), nrow = P, ncol = V)))
  D <- (1 - G) * s2 * E + G * t *l2
  b <- fast_sampler(Phi = X / sqrt(s2),
                    D   = D,
                    a   = y / sqrt(s2))
  
  #############################################################################
  # Horse Shoe Structure
  #############################################################################
  # Samples l2
  l2 <- 1 / rgamma(n     = p,
                   shape = 1,
                   rate  = 1 / v + (G * b^2) / (2 * t * s2))
  
  # Samples t2
  for(i in 1:P){
    if(g[i] == 1){
      t2[i] <- 1 / rgamma(n     = 1,
                          shape = (p + 1) / 2,
                          rate  = 1 / xi + (1 / (2 * s2)) * sum(b^2 / l2))
    } else {
      t2[i] <- 1 / rgamma(n     = 1,
                          shape = 1,
                          rate  = 1 / xi)
    }
  }
  
  # Samples v2
  v  <- 1 / rgamma(n     = p,
                   shape = 1,
                   rate  = 1 + 1 / l2)
  
  # Samples xi
  xi <- 1 / rgamma(n     = 1,
                   shape = 1,
                   rate  = 1 + 1 / t2)
  
  # Samples e2
  B <- matrix(data = B, nrow = V, ncol = P)
  for(i in 1:P){
    if(g[i] == 0){
      e2[i] <- 1 / rgamma(n     = 1,
                          shape = V/2 + 1,
                          rate  = crossprod(B[, i])/ (2 * s2) + 1)
    } else {
      e2[i] <- 1 / rgamma(n     = 1,
                          shape = 1,
                          rate  = 0.01)
    }
  }
  
  # Samples g
  L <- matrix(data = l2, nrow = V, ncol = P)
  for(i in 1:P){
    print(i)
    p0   <- dnorm(x    = B[,i],
                  mean = 0,
                  sd   = sqrt(s2 * e2),
                  log  = TRUE)
    print(p0)
    p1   <- dnorm(x    = B[,i],
                  mean = 0,
                  sd   = sqrt(s2 * t2[i] * L[, i]),
                  log  = TRUE)
    print(p1)
    p1   <- exp(sum(p1) - sum(p0))
    if(is.finite(p1)){
      p1 <- p1 / (p1 + 1)
    } else {
      p1 <- 1
    }
    print(p1)
    g[i] <- rbinom(n    = 1,
                   size = 1,
                   prob = p1)
  }
  
  #############################################################################
  # Samples s2
  #############################################################################
  G <- c(t(matrix(data = rep(g, V), nrow = P, ncol = V)))
  t <- c(t(matrix(data = rep(t2, V), nrow = P, ncol = V)))
  D <- (1 - G) * s2 * e2 + G * t *l2
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + p) / 2,
                   rate  = crossprod(x = y - X %*% b) / 2 +
                     t(b / D) %*% b / 2)
  
  #############################################################################
  # Returns the Samples
  #############################################################################
  return(list(b  = b,
              s2 = s2,
              l2 = l2,
              t2 = t2,
              v  = v,
              xi = xi,
              e2 = e2,
              g  = g))
}