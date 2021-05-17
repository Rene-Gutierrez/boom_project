# Fast Horse Shoe Sampler (For Big p and small n)

su_boom_iteration <- function(y,
                              X,
                              t2, # Horseshoe Structure for B
                              l2, #
                              xi, #
                              v,  #
                              c2, # Spike & Slab Structure for Theta
                              m2, #
                              gg, #
                              s2,
                              e2){
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
  # Generates g
  g     <- 1 - apply(X = 1 - gg, MARGIN = 1, FUN = prod)
  # Prior Diagonal
  temp1 <- s2 * c(t2 * t(t(l2) * g + (1 - g) * e2))
  temp2 <- (c2 * m2 * gg + (1 - gg) * m2)[lower.tri(gg)]
  Lam   <- c(temp1, temp2)
  # Fast Sampling
  b <- fast_sampler(Phi = X / sqrt(s2),
                    D   = Lam,
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
  # Spike & Slab Structure for Theta
  #############################################################################
  
  #############################################################################
  # Samples s2
  #############################################################################
  temp1 <- c(t2 * t(t(l2) * g + (1 - g) * e2))
  Lam   <- c(temp1)
  s2 <- 1 / rgamma(n     = 1,
                   shape = (n + P * V) / 2,
                   rate  = crossprod(x = y - X %*% b) / 2 +
                     t(c(B)[Lam != 0] / Lam[Lam != 0]) %*% c(B)[Lam != 0] / 2)
  
  #############################################################################
  # Samples g
  #############################################################################
  pg <- rep(0, P)
  for(i in 1:P){
    # print(i)
    # print(g)
    # print('p0')
    # print(sum(dnorm(x    = B[, i],
    #             mean = 0,
    #             sd   = sqrt(s2 * e2),
    #             log  = TRUE)))
    # print(sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
    #                 mean = 0,
    #                 sd   = sqrt(c2 * m2),
    #                 log  = TRUE)))
    # print(dnorm(x    = Theta[i, -i][g[-i] == 1],
    #             mean = 0,
    #             sd   = sqrt(s2 * o2),
    #             log  = TRUE))
    # print(sum(dnorm(x    = B[, i],
    #                 mean = 0,
    #                 sd   = sqrt(s2 * e2),
    #                 log  = TRUE)) + log(1/2) +
    #         sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
    #                   mean = 0,
    #                   sd   = sqrt(s2 * o2),
    #                   log  = TRUE)))
    p0   <- sum(dnorm(x    = B[, i],
                      mean = 0,
                      sd   = sqrt(s2 * e2),
                      log  = TRUE)) + log(1/2) +
      sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
                mean = 0,
                sd   = sqrt(m2),
                log  = TRUE))
    # print('p1')
    # print(sum(dnorm(x    = B[, i],
    #             mean = 0,
    #             sd   = sqrt(s2 * t2 * l2[, i]),
    #             log  = TRUE)))
    # print(sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
    #                 mean = 0,
    #                 sd   = sqrt(s2 * r2 * m2[, i]),
    #                 log  = TRUE)))
    # print(dnorm(x    = Theta[i, -i][g[-i] == 1],
    #             mean = 0,
    #             sd   = sqrt(s2 * r2 * m2[-i, i][g[-i] == 1]),
    #             log  = TRUE))
    # print(sum(dnorm(x    = B[, i],
    #                 mean = 0,
    #                 sd   = sqrt(s2 * t2 * l2[, i]),
    #                 log  = TRUE)) + log(1/2) +
    #         sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
    #                   mean = 0,
    #                   sd   = sqrt(s2 * r2 * m2[, i]),
    #                   log  = TRUE)))
    p1   <- sum(dnorm(x    = B[, i],
                      mean = 0,
                      sd   = sqrt(s2 * t2 * l2[, i]),
                      log  = TRUE)) + log(1/2) +
      sum(dnorm(x    = Theta[i, -i][g[-i] == 1],
                mean = 0,
                sd   = sqrt(c2 * m2),
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
  # g <- c(rep(1, QQ), rep(0, P - QQ))
  # g <- rep(1, P)
  
  temp1 <- s2 * c(t2 * t(t(l2) * g + (1 - g) * e2))
  temp2 <- (c2 * m2 * (g %*% t(g)) + m2 * (1 - (g %*% t(g))))[lower.tri((g %*% t(g)))]
  Lam   <- c(temp1, temp2)
  #############################################################################
  # Returns the Samples
  #############################################################################
  return(list(B      = B,
              l2     = l2,
              v      = v,
              t2     = t2,
              xi     = xi,
              Theta  = Theta,
              s2     = s2,
              g      = g,
              pg     = pg,
              D      = Lam,
              flag   = flag))
}