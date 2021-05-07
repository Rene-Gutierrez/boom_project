# Fast Horse Shoe Sampler (For Big p and small n)

u_boom_iteration <- function(y,
                             X,
                             t2, # Horseshoe Structure for B
                             l2, #
                             xi, #
                             v,  #
                             r2, # Horseshoe Structure for Theta
                             m2, #
                             nu, #
                             w,  #
                             s2,
                             g,
                             e2,
                             o2){
  #############################################################################
  # Set-Up
  #############################################################################
  n <- dim(X)[1]
  V <- dim(l2)[1]
  P <- dim(l2)[2]
  
  # Numerical Failure Flag
  flag <- FALSE
  
  #############################################################################
  # Samples B and Theta
  #############################################################################
  # Prior Diagonal
  Lam <- c(t2 * c(l2), r2 * c(m2[lower.tri(m2)]))
  # Fast Sampling
  b <- fast_sampler(Phi = X / sqrt(s2),
                    D   = s2 * Lam,
                    a   = y / sqrt(s2))
  # Distributes the Estimate
  B     <- matrix(data = b[1:(V * P)], nrow = V, ncol = P)
  Theta <- diag(P)
  Theta[upper.tri(Theta, diag = TRUE)] <- 0
  Theta[lower.tri(Theta)] <- b[(V * P + 1):(V * P + P * (P - 1) / 2)]
  Theta <- Theta + t(Theta)
  #############################################################################
  # Horse Shoe Structure for B
  #############################################################################
  # Samples l2
  for(i in 1:P){
    if(g[i] == 1){
      l2[, i] <- 1 / rgamma(n     = V,
                            shape = 1,
                            rate  = 1 / v[, i] + B[, i]^2 / (2 * s2 * t2))
    } else {
      l2[, i] <- 1 / rgamma(n     = V,
                            shape = 1 / 2,
                            rate  = 1 / v[, i])
    }
  }
  
  # Samples t2
  if(sum(g) > 0){
    t2 <- 1 / rgamma(n     = 1,
                     shape = (sum(g) * V + 1) / 2,
                     rate  = 1 / xi +
                       sum(Lambda(M = B^2, g = g, a = 0) / (2 * s2 * l2)))
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
  
  
  #############################################################################
  # Horse Shoe Structure for Theta
  #############################################################################
  # Samples m2
  for(i in 2:P){
    for(j in 1:(i - 1)){
      if(g[i] * g[j] == 1){
        m2[i, j] <- 1 / rgamma(n     = 1,
                               shape = 1,
                               rate  = 1 / w[i, j] + Theta[i, j]^2 / (2 * s2 * r2))
        m2[j, i] <- m2[i, j] 
      } else {
        m2[i, j] <- 1 / rgamma(n     = 1,
                               shape = 1 / 2,
                               rate  = 1 / w[i, j])
        m2[j, i] <- m2[i, j]
      }
    }
  }
  
  # Samples r2
  temp1 <- Theta[g == 1, g == 1]
  temp1 <- temp1[lower.tri(temp1)]
  temp2 <- m2[g == 1, g == 1]
  temp2 <- temp2[lower.tri(temp2)]
  if(sum(g) > 0){
    r2 <- 1 / rgamma(n     = 1,
                     shape = (sum(g) * (sum(g) - 1) / 2 + 1) / 2,
                     rate  = 1 / nu +
                       sum(temp1^2 / (2 * s2 * temp2)))
  } else {
    r2 <- 1 / rgamma(n     = 1,
                     shape = 1 / 2,
                     rate  = 1 / nu)
  }
  
  # Samples w
  out <- 1 / rgamma(n     = P * (P - 1) / 2,
                    shape = 1,
                    rate  = 1 + 1 / m2[lower.tri(m2)])
  diag(w)         <- 0
  w[lower.tri(w)] <- out
  w               <- w + t(w)
  
  # Samples nu
  nu <- 1 / rgamma(n     = 1,
                   shape = 1,
                   rate  = 1 + 1 / r2)
  
  #############################################################################
  # Samples s2
  #############################################################################
  Lam <- c(t2 * c(l2), r2 * c(m2[lower.tri(m2)]))
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + P * V + P * (P - 1) / 2) / 2,
                   rate  = crossprod(x = y - X %*% b) / 2 + t(b / Lam) %*% b / 2)
  
  #############################################################################
  # Samples g
  #############################################################################
  pg <- rep(0, P)
  # for(i in 1:P){
  #   p0   <- sum(dnorm(x    = B[, i],
  #                     mean = 0,
  #                     sd   = sqrt(s2 * e2),
  #                     log  = TRUE)) + log(1/2) + 
  #     sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
  #               mean = 0,
  #               sd   = m2,
  #               log  = TRUE))
  #   p1   <- sum(dnorm(x    = B[, i],
  #                     mean = 0,
  #                     sd   = sqrt(s2 * t2 * l2[, i]),
  #                     log  = TRUE)) + log(1/2) + 
  #     sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
  #               mean = 0,
  #               sd   = M2,
  #               log  = TRUE))
  #   p1   <- exp(sum(p1) - sum(p0))
  #   if(is.finite(p1)){
  #     p1 <- p1 / (p1 + 1)
  #   } else {
  #     p1 <- 1
  #   }
  #   g[i] <- rbinom(n    = 1,
  #                  size = 1,
  #                  prob = p1)
  #   pg[i] <- p1
  # }
  g <- rep(1, P)
  
  #############################################################################
  # Returns the Samples
  #############################################################################
  return(list(B      = B,
              l2     = l2,
              v      = v,
              t2     = t2,
              xi     = xi,
              Theta  = Theta,
              m2     = m2,
              w      = w,
              r2     = r2,
              nu     = nu,
              s2     = s2,
              g      = g,
              pg     = pg,
              flag   = flag))
}