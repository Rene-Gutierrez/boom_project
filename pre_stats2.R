pre_stats2 <- function(ty,
                      A,
                      G,
                      B,
                      Theta){
  
  # Dimensions
  n <- length(ty)
  S <- dim(B)[1]
  V <- dim(B)[2]
  P <- dim(B)[3]
  
  # Computes the Predictive Response
  AX      <- A
  dim(AX) <- c(n, P * P)
  AX      <- AX[, lower.tri(Theta[1,,])]
  GX      <- G
  dim(GX) <- c(n, V * P)
  vB      <- B
  dim(vB) <- c(S, V * P)
  vT      <- Theta
  dim(vT) <- c(S, P * P)
  vT      <- vT[, lower.tri((Theta[1,,]))]
  X       <- cbind(GX, AX)
  b       <- cbind(vB, vT)
  y       <- X %*% t(b)
  
  # Coverage
  low <- apply(y, 1, quantile, pro = 0.025)
  upp <- apply(y, 1, quantile, pro = 0.975)
  cov <- mean((low <= ty) * (upp >= ty))
  
  # Length
  len <- mean(upp - low)
  
  # Mse
  vB      <- apply(X = vB, MARGIN = 2, FUN = median)
  vT      <- Theta
  dim(vT) <- c(S, P * P)
  vT      <- vT[, lower.tri((Theta[1,,]))]
  vT      <- apply(X = vT, MARGIN = 2, FUN = median)
  X       <- cbind(GX, AX)
  b       <- c(vB, vT)
  y       <- X %*% b
  mse  <- mean((y - ty)^2)
  pmse <- mse / var(ty)
  
  return(list(y    = y,
              err  = colMeans((y - ty)^2),
              mse  = mse,
              pmse = pmse,
              cov  = cov,
              len  = len))
}