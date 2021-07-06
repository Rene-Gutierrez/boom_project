tra_pre_stats <- function(ty,
                          A,
                          G,
                          coe){
  
  # Dimesions
  n <- length(ty)
  V <- dim(G)[2]
  P <- dim(G)[3]
  
  # Computes the Predictive Response
  AX      <- A
  dim(AX) <- c(n, P * P)
  AX      <- AX[, lower.tri(A[1,,])]
  GX      <- G
  dim(GX) <- c(n, V * P)
  X       <- cbind(GX, AX)
  b       <- coe
  y       <- X %*% b
  
  # Mse
  mse  <- mean((y - ty)^2)
  pmse <- mse / var(ty) 
  
  # Coverage
  cov <- NA
  
  # Length
  len <- NA
  
  return(list(mse  = mse,
              pmse = pmse,
              cov  = cov,
              len  = len))
}