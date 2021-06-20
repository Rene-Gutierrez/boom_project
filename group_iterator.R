group_iterator <- function(AX,
                           GX,
                           y,
                           Theta,
                           LT,
                           B,
                           LB,
                           s2,
                           g,
                           r,
                           p,
                           tM,
                           bM){
  # Problem dimensions
  n <- length(y)  # Number of Observations
  V <- dim(B)[1]  # Voxel Size
  P <- dim(B)[2]  # Number of ROI's
  q <- sum(g[-p]) # Number of Active Regions (Excluding the Current Region)
  
  # Auxiliary Variables
  ipb <- bM[, p]
  ipt <- tM[p, -p][g[-p] == 1]
  iqb <- c(bM[, -p][, g[-p] == 1])
  iqt <- tM[-p, -p][g[-p] == 1, g[-p] == 1]
  iqt <- iqt[lower.tri(iqt)]
  
  # Cases depending on the number of Active regions
  if(q == 0){
    # Auxiliary Covariate
    X <- GX[,ipb]
    # Auxiliary Response
    R <- y
    # Obtains the Point Estimates
    LL <- diag(1 / LB[,p])
    Q  <- t(X) %*% X + LL 
    U  <- chol(Q)
    bh <- backsolve(U, forwardsolve(t(U), t(X) %*% R))
    # Odds for g == 1
    od <- - log(det(Q)) / 2
    od <- od + t(bh) %*% Q %*% bh / s2 / 2
    od <- od - log(prod(diag(LL))) / 2
    od <- od + log(r) - log(1 - r)
    od <- exp(od)
    # Probability for g == 1
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    # Updates g
    g[p] <- rbinom(n = 1, size = 1, prob = pro)
    # Updates B
    if(g[p] == 1){
      b      <- backsolve(U, rnorm(n = V)) + bh
      B[, p] <- b
    } else {
      b      <- rep(0, V)
      B[, p] <- b
    }
  } else if(q == 1){
    # Auxiliary Covariate
    X1 <- GX[, ipb]
    X2 <- AX[, ipt]
    X  <- cbind(X1, X2)
    # Auxiliary Response
    Z <- GX[, iqb]
    R <- y - Z %*% c(B[,-p][, g[-p] == 1])
    # Obtains the Point Estimates
    LL <- diag(1 / c(LB[, p], LT[p,-p][g[-p] == 1]))
    Q  <- t(X) %*% X + LL 
    U  <- chol(Q)
    bh <- backsolve(U, forwardsolve(t(U), t(X) %*% R))
    # Odds for g == 1
    od <- - log(det(Q)) / 2
    od <- od + t(bh) %*% Q %*% bh / s2 / 2
    od <- od - log(prod(diag(LL))) / 2
    od <- od + log(r) - log(1 - r)
    od <- exp(od)
    # Probability for g == 1
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    # Updates g
    g[p] <- rbinom(n = 1, size = 1, prob = pro)
    # Updates B and Theta
    if(g[p] == 1){
      b      <- backsolve(U, rnorm(n = V + q)) + bh
      B[, p] <- b[1:V]
      Theta[p, -p][g[-p] == 1] <- b[(V + 1):(V + q)]
      Theta[-p, p][g[-p] == 1] <- b[(V + 1):(V + q)]
    } else {
      b      <- rep(0, V + q)
      B[, p] <- b[1:V]
      Theta[p, -p][g[-p] == 1] <- b[(V + 1):(V + q)]
      Theta[-p, p][g[-p] == 1] <- b[(V + 1):(V + q)]
    }
  } else {
    # Auxiliary Covariate
    X1 <- GX[, ipb]
    X2 <- AX[, ipt]
    X  <- cbind(X1, X2)
    # Auxiliary Response
    Z1 <- GX[, iqb]
    Z2 <- AX[, iqt]
    Z  <- cbind(Z1, Z2)
    R  <- y - Z %*% c(c(B)[iqb], Theta[lower.tri(Theta)][iqt])
    # Obtains the Point Estimates
    LL <- diag(1 / c(LB[, p], LT[p,-p][g[-p] == 1]))
    Q  <- t(X) %*% X + LL
    U  <- chol(Q)
    bh <- backsolve(U, forwardsolve(t(U), t(X) %*% R))
    # Odds for g == 1
    od <- - log(det(Q)) / 2
    od <- od + t(bh) %*% Q %*% bh / s2 / 2
    od <- od - log(prod(diag(LL))) / 2
    od <- od + log(r) - log(1 - r)
    od <- exp(od)
    # Probability for g == 1
    if(is.infinite(od)){
      pro <- 1
    } else {
      pro <- od / (1 + od)
    }
    # Updates g
    g[p] <- rbinom(n = 1, size = 1, prob = pro)
    # Updates B and Theta
    if(g[p] == 1){
      b      <- backsolve(U, rnorm(n = V + q)) + bh
      B[, p] <- b[1:V]
      Theta[p, -p][g[-p] == 1] <- b[(V + 1):(V + q)]
      Theta[-p, p][g[-p] == 1] <- b[(V + 1):(V + q)]
    } else {
      b      <- rep(0, V + q)
      B[, p] <- b[1:V]
      Theta[p, -p][g[-p] == 1] <- b[(V + 1):(V + q)]
      Theta[-p, p][g[-p] == 1] <- b[(V + 1):(V + q)]
    }
  }
  return(list(bh    = bh,
              b     = b,
              od    = od,
              pro   = pro,
              ipb   = ipb,
              ipt   = ipt,
              iqb   = iqb,
              iqt   = iqt,
              R     = R,
              LL    = LL,
              U     = U,
              B     = B,
              Theta = Theta,
              g     = g))
}