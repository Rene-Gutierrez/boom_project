tra_sta <- function(tTheta,
                    tB,
                    tg,
                    tgB,
                    coe){
  # Dimensions
  V <- dim(tB)[1]
  P <- dim(tB)[2]
  
  # Recreates Theta or B
  B     <- matrix(data = coe[1:(V * P)], nrow = V, ncol = P)
  Theta <- coe[(V * P + 1):(V * P + P * (P - 1) / 2)]
  tem   <- matrix(data = 0, nrow = P, ncol = P)
  tem[lower.tri(tem)] <- Theta
  tem <- tem + t(tem)
  
  # Creates g
  selB <- colSums(B != 0) > 0
  selT <- colSums(tem != 0) > 0
  g    <- (selB * selT) + 0 
  
  # TPR and TNR
  TPR <- mean(g[tg == 1])
  TNR <- 1 - mean(g[tg == 0])
  
  # MSE Non-Zero Coefficients Theta
  vT      <- Theta
  vtT     <- tTheta[lower.tri(tTheta)]
  selT    <- ((tg %*% t(tg)) == 1)
  selT    <- selT[lower.tri(selT)]
  MseNzT  <- mean((t(vT[selT]) - vtT[selT])^2) 
  # MSE Zero Coefficients Theta
  MsezT  <- mean((t(vT[!selT]) - vtT[!selT])^2)
  # MSE Theta
  MseT <- mean((vT - vtT)^2)
  
  # MSE Global and Local Non-Zero Coefficients B
  vB      <- B
  dim(vB) <- c(V * P)
  MseNzB  <- mean((t(vB[c(tgB) == 1]) - tB[tgB == 1])^2)
  # MSE Local Zero
  vNzB      <- B[,tg == 1]
  dim(vNzB) <- c(V * sum(tg == 1)) 
  MselzB    <- mean((vNzB[tgB[, tg == 1] == 0])^2)
  # MSE Global Zero
  vNzB      <- B[,tg == 0] 
  MsegzB    <- mean(vNzB^2)
  # MSE B
  MseB <- mean((vB - c(tB))^2)
  
  # Global MSE
  Mse <- mean(c((vT - vtT)^2, (vB - c(tB))^2))
  
  return(list(B     = B,
              Theta = Theta,
              selB  = selB,
              selT  = selT,
              g     = g,
              TPR   = TPR,
              TNR   = TNR,
              MseNzT = MseNzT,
              MsezT  = MsezT,
              MseT   = MseT,
              MseNzB = MseNzB,
              MselzB = MselzB,
              MsegzB = MsegzB,
              MseB   = MseB,
              Mse    = Mse))
}