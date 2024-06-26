### Spike & Slab Point Mass Structure

groupSymSpikeSlabUpd <- function(A, y, Theta, L, s, g, p, r){
  n      <- length(y)
  P      <- dim(A)[2]
  if(sum(g[-p]) == 0){
    print('Case g-p = 0')
    od  <- r / (1 - r)
    pro <-  od / (1 + od)
    ###########################################################################
    # Samples g
    ###########################################################################
    g[p]  <- rbinom(n = 1, size = 1, prob = pro)
    Th    <- c()
    Ch    <- c()
  } else if(sum(g[-p]) == 1) {
    print('Case g-p = 1')
    ###########################################################################
    # Covariate
    ###########################################################################
    X      <- A[, p, -p][, g[-p] == 1]
    ###########################################################################
    # Obtains R
    ###########################################################################
    R      <- y
    ###########################################################################
    # Obtains Theta Hat
    ###########################################################################
    LL <- L[p, -p][g[-p] == 1]
    LL <- 1 / LL
    Th <- t(X) %*% R / (t(X) %*% X + LL)
    ###########################################################################
    # Computes the Odds for g = 1
    ###########################################################################
    od <- - log(t(X) %*% X + LL) / 2
    od <- od + Th * (t(X) %*% X + LL) * Th / s / 2
    od <- od - log(LL) / 2
    od <- od + log(r) - log(1 - r)
    od <- exp(od)
    ###########################################################################
    # Computes Probability for g = 1
    ###########################################################################
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    ###########################################################################
    # Updates g and Theta
    ###########################################################################
    g[p] <- rbinom(n = 1, size = 1, prob = pro)
    if(g[p] == 1){
      Ch <- sqrt(t(X) %*% X + LL)
      Tu <- rnorm(n = 1, mean = Th, sd = 1 / sqrt(t(X) %*% X + LL))
      Theta[p, -p][g[-p] == 1] <- Tu
      Theta[-p, p][g[-p] == 1] <- Tu
    } else {
      Ch <- matrix(data = NA, nrow = sum(g[-p]), ncol = sum(g[-p]))
      Theta[p, -p][g[-p] == 1] <- 0
      Theta[-p, p][g[-p] == 1] <- 0
    }
  } else if(sum(g[-p]) == 2){
    print('Case g-p = 2')
    ###########################################################################
    # Covariate
    ###########################################################################
    X      <- A[, p, -p][, g[-p] == 1]
    ###########################################################################
    # Obtains R
    ###########################################################################
    Z      <- A[, -p, -p][, g[-p] == 1, g[-p] == 1]
    dim(Z) <- c(n, sum(g[-p]) * sum(g[-p]))
    Z      <- Z[, lower.tri(diag(sum(g[-p])))]
    W      <- Theta[-p, -p][g[-p] == 1, g[-p] == 1][lower.tri(diag(sum(g[-p])))]
    R      <- y - Z * W
    ###########################################################################
    # Obtains Theta Hat
    ###########################################################################
    LL <- L[p, -p][g[-p] == 1]
    LL <- diag(1 / LL)
    Th <- solve(t(X) %*% X + LL, t(X) %*% R)
    ###########################################################################
    # Computes the Odds for g = 1
    ###########################################################################
    od <- - log(det(t(X) %*% X + LL)) / 2
    od <- od + t(Th) %*% (t(X) %*% X + LL) %*% Th / s / 2
    od <- od - log(prod(diag(LL))) / 2
    od <- od + log(r) - log(1 - r)
    od <- exp(od)
    ###########################################################################
    # Computes Probability for g = 1
    ###########################################################################
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    ###########################################################################
    # Updates g and Theta
    ###########################################################################
    g[p] <- rbinom(n = 1, size = 1, prob = pro)
    if(g[p] == 1){
      Ch <- chol((t(X) %*% X + LL))
      Tu <- solve(t(Ch), rnorm(n = 2)) + Th
      Theta[p, -p][g[-p] == 1] <- Tu
      Theta[-p, p][g[-p] == 1] <- Tu
    } else {
      Ch <- matrix(data = NA, nrow = sum(g[-p]), ncol = sum(g[-p]))
      Theta[p, -p][g[-p] == 1] <- 0
      Theta[-p, p][g[-p] == 1] <- 0
    }
  } else {
    print('General Case')
    ###########################################################################
    # Covariate
    ###########################################################################
    X      <- A[, p, -p][, g[-p] == 1]
    ###########################################################################
    # Obtains R
    ###########################################################################
    Z      <- A[, -p, -p][, g[-p] == 1, g[-p] == 1]
    dim(Z) <- c(n, sum(g[-p]) * sum(g[-p]))
    Z      <- Z[, lower.tri(diag(sum(g[-p])))]
    W      <- Theta[-p, -p][g[-p] == 1, g[-p] == 1][lower.tri(diag(sum(g[-p])))]
    R      <- y - Z %*% W
    ###########################################################################
    # Obtains Theta Hat
    ###########################################################################
    LL <- L[p, -p][g[-p] == 1]
    LL <- diag(1 / LL)
    Th <- solve(t(X) %*% X + LL, t(X) %*% R)
    ###########################################################################
    # Computes the Odds for g = 1
    ###########################################################################
    od <- - log(det(t(X) %*% X + LL)) / 2
    od <- od + t(Th) %*% (t(X) %*% X + LL) %*% Th / s / 2
    od <- od - log(prod(diag(LL))) / 2
    od <- od + log(r) - log(1 - r)
    od <- exp(od)
    ###########################################################################
    # Computes Probability for g = 1
    ###########################################################################
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    ###########################################################################
    # Updates g and Theta
    ###########################################################################
    g[p] <- rbinom(n = 1, size = 1, prob = pro)
    if(g[p] == 1){
      Ch <- chol((t(X) %*% X + LL))
      Tu <- solve(t(Ch), rnorm(n = sum(g[-p]))) + Th
      Theta[p, -p][g[-p] == 1] <- Tu
      Theta[-p, p][g[-p] == 1] <- Tu
    } else {
      Ch <- matrix(data = NA, nrow = sum(g[-p]), ncol = sum(g[-p]))
      Theta[p, -p][g[-p] == 1] <- 0
      Theta[-p, p][g[-p] == 1] <- 0
    }
  }
  return(list(pro   = pro,
              g     = g,
              Theta = Theta,
              Th    = Th,
              Ch    = Ch))
}