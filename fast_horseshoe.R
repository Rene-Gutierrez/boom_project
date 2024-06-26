# Fast Horse Shoe Sampler (For Big p and small n)

fast_horseshoe <- function(y,
                           X,
                           xi,
                           v,
                           l2,
                           t2,
                           s2){
  #############################################################################
  # Dimensions
  #############################################################################
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #############################################################################
  # Samples b (Fast)
  #############################################################################
  b <- fast_sampler(Phi = X / sqrt(s2),
                    D   = s2 * t2 * l2,
                    a   = y / sqrt(s2))
  
  #############################################################################
  # Horse Shoe Structure
  #############################################################################
  # Samples l2
  l2 <- 1 / rgamma(n     = p,
                   shape = 1,
                   rate  = 1 / v + b^2 / (2 * t2 * s2))
  
  # Samples t2
  t2 <- 1 / rgamma(n     = 1,
                   shape = (p + 1) / 2,
                   rate  = 1 / xi + (1 / (2 * s2)) * sum(b^2 / l2))
  
  # Samples v2
  v  <- 1 / rgamma(n     = p,
                   shape = 1,
                   rate  = 1 + 1 / l2)
  
  # Samples xi
  xi <- 1 / rgamma(n     = 1,
                   shape = 1,
                   rate  = 1 + 1 / t2)
  
  #############################################################################
  # Samples s2
  #############################################################################
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + p) / 2,
                   rate  = crossprod(x = y - X %*% b) / 2 +
                     t(b / (t2 * l2)) %*% b / 2)
  
  #############################################################################
  # Returns the Samples
  #############################################################################
  return(list(b  = b,
              s2 = s2,
              l2 = l2,
              t2 = t2,
              v  = v,
              xi = xi))
}