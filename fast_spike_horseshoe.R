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
  # Set-Up
  #############################################################################
  n <- dim(X)[1]
  V <- dim(l2)[1]
  P <- dim(l2)[2]
  
  #############################################################################
  # Samples b (Fast)
  #############################################################################
  b <- fast_sampler(Phi = X / sqrt(s2),
                    D   = c(t((1 - g) * e2 + t(l2) * g * t2)),
                    a   = y / sqrt(s2))
  b <- matrix(data = b, nrow = V, ncol = P)
  #############################################################################
  # Horse Shoe Structure
  #############################################################################
  # Samples l2
  l2 <- 1 / rgamma(n     = P * V,
                   shape = 1,
                   rate  = 1 / c(v) + c(t(t(b^2) * g / (2 * t2 * s2))))
  l2 <- matrix(data = l2, nrow = V, ncol = P)
  
  # Samples t2
  for(i in 1:P){
    if(g[i] == 1){
      t2[i] <- 1 / rgamma(n     = 1,
                          shape = (V + 1) / 2,
                          rate  = 1 / xi[i] +
                            (1 / (2 * s2)) * sum(b[, i]^2 / l2[, i]))
    } else {
      t2[i] <- 1 / rgamma(n     = 1,
                          shape = 1,
                          rate  = 1 / xi[i])
    }
  }
  
  # Samples v2
  v  <- 1 / rgamma(n     = P * V,
                   shape = 1,
                   rate  = 1 + 1 / l2)
  v <- matrix(data = v, nrow = V, ncol = P)
  
  # Samples xi
  xi <- 1 / rgamma(n     = P,
                   shape = 1,
                   rate  = 1 + 1 / t2)
  
  # Samples e2
  for(i in 1:P){
    if(g[i] == 0){
      e2[i] <- 1 / rgamma(n     = 1,
                          shape = V/2 + 1,
                          rate  = crossprod(b[, i])/ (2 * s2) + 1)
    } else {
      e2[i] <- 1 / rgamma(n     = 1,
                          shape = 1,
                          rate  = 0.01)
    }
  }
  e2 <- rep(0.01, P)
  
  # Samples g
  pg <- rep(0, P)
  for(i in 1:P){
    p0   <- dnorm(x    = b[,i],
                  mean = 0,
                  sd   = sqrt(s2 * e2),
                  log  = TRUE)
    p1   <- dnorm(x    = b[,i],
                  mean = 0,
                  sd   = sqrt(s2 * t2[i] * l2[, i]),
                  log  = TRUE)
    p1   <- exp(sum(p1) - sum(p0))
    if(is.finite(p1)){
      p1 <- p1 / (p1 + 1)
    } else {
      p1 <- 1
    }
    g[i] <- rbinom(n    = 1,
                   size = 1,
                   prob = p1)
    pg[i] <- p1
  }
  
  #############################################################################
  # Samples s2
  #############################################################################
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + P * V) / 2,
                   rate  = crossprod(x = y - X %*% c(b)) / 2 +
                     t(c(b) / c(t((1 - g) * e2 + t(l2) * g * t2))) %*% c(b) / 2)
  
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
              g  = g,
              pg = pg))
}