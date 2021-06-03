# Fast Horse Shoe Sampler (For Big p and small n)

u_boom_iteration <- function(y,
                             X,
                             xi,
                             v,
                             l2,
                             t2,
                             s2,
                             g,
                             e2,
                             h2,
                             m2,
                             M2,
                             K){
  #############################################################################
  # Set-Up
  #############################################################################
  n <- dim(X)[1]
  V <- dim(l2)[1]
  P <- dim(l2)[2]
  
  # Numerical Failure Flag
  flag <- FALSE
  
  #############################################################################
  # Samples b (Fast)
  #############################################################################
  out <- fast_sampler(Phi = X / sqrt(s2),
                      D   = L(M = s2 * t2 * l2, g = g, a = s2 * e2, m2 = m2, M2 = M2),
                      a   = y / sqrt(s2)) 
  b   <- matrix(data = out[1:(V * P)], nrow = V, ncol = P)
  Theta <- diag(P)
  Theta[upper.tri(Theta, diag = TRUE)] <- 0
  Theta[lower.tri(Theta)] <- out[(V * P + 1):(V * P + P * (P - 1) / 2)]
  Theta <- Theta + t(Theta)
  #############################################################################
  # Horse Shoe Structure
  #############################################################################
  # Samples l2
  for(i in 1:P){
    if(g[i] == 1){
      l2[, i] <- 1 / rgamma(n     = V,
                            shape = 1,
                            rate  = 1 / v[,i] + b[,i]^2 / (2 * s2 * t2))
    } else {
      l2[, i] <- 1 / rgamma(n     = V,
                            shape = 1 / 2,
                            rate  = 1 / v[,i])
    }
  }
  
  # # Samples e2
  if(sum(1 - g) > 0){
    e2 <- 1 / rgamma(n     = 1,
                     shape = (sum(1 - g) * V + 1) / 2,
                     rate  = 1 / h2 + 
                       sum(Lambda(M = b^2, g = 1 - g, a = 0) / (2 * s2)))
  } else {
    e2 <- 1 / rgamma(n     = 1,
                     shape = 1 / 2,
                     rate  = 1 / h2)
  }
  
  # Samples t2
  if(sum(g) > 0){
    t2 <- 1 / rgamma(n     = 1,
                     shape = (sum(g) * V + 1) / 2,
                     rate  = 1 / xi +
                       sum(Lambda(M = b^2, g = g, a = 0) / (2 * s2 * l2)))
  } else {
    t2 <- 1 / rgamma(n     = 1,
                     shape = 1 / 2,
                     rate  = 1 / xi)
  }
  
  # Samples v2
  v  <- 1 / rgamma(n     = P * V,
                   shape = 1,
                   rate  = 1 + 1 / c(l2))
  v <- matrix(data = v, nrow = V, ncol = P)
  
  # Samples xi
  xi <- 1 / rgamma(n     = 1,
                   shape = 1,
                   rate  = 1 + 1 / t2)
  
  # Samples h2
  h2 <- 1 / rgamma(n     = 1,
                   shape = 1,
                   rate  = 1 + 1 / e2)
  
  #############################################################################
  # Samples Theta
  #############################################################################
  
  
  #############################################################################
  # Samples s2
  #############################################################################
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + P * V) / 2,
                   rate  = crossprod(x = y - X %*% c(c(b), c(Theta[lower.tri(Theta)]))) / 2 +
                     t(c(b) / c(Lambda(M = l2 * t2, g = g, a = e2))) %*% c(b) / 2)
  
  
  #############################################################################
  # Samples g
  #############################################################################
  pg <- rep(0, P)
  for(i in 1:P){
    p0   <- sum(dnorm(x    = b[, i],
                      mean = 0,
                      sd   = sqrt(s2 * e2),
                      log  = TRUE)) + log(1/2) + 
            sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
                      mean = 0,
                      sd   = m2,
                      log  = TRUE))
    p1   <- sum(dnorm(x    = b[, i],
                      mean = 0,
                      sd   = sqrt(s2 * t2 * l2[, i]),
                      log  = TRUE)) + log(1/2) + 
            sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
                      mean = 0,
                      sd   = M2,
                      log  = TRUE))
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
  # g <- rep(1, P)
  
  #############################################################################
  # Returns the Samples
  #############################################################################
  return(list(b      = b,
              s2     = s2,
              l2     = l2,
              t2     = t2,
              v      = v,
              xi     = xi,
              e2     = e2,
              h2     = h2,
              g      = g,
              pg     = pg,
              Theta  = Theta,
              flag   = flag))
}