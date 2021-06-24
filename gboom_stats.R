gboom_stats <- function(tTheta,
                        tB,
                        tg,
                        tgB,
                        Theta,
                        B,
                        g      = NA,
                        flag   = TRUE){
  # Dimensions
  S <- dim(B)[1]
  V <- dim(B)[2]
  P <- dim(B)[3]
  
  if(flag){
    # True Positive Rate
    TPR <- mean(g[, tg])
    TNR <- 1 - mean(g[, tg == 0])
  } else {
    TPR <- NA
    TNR <- NA
  }
  
  
  # MSE Non-Zero Coefficients Theta
  vT      <- Theta
  dim(vT) <- c(S, P * P)
  vT      <- vT[, lower.tri(Theta[1,,])]
  vtT     <- tTheta[lower.tri(tTheta)]
  selT    <- ((tg %*% t(tg)) == 1)
  selT    <- selT[lower.tri(selT)]
  MseNzT  <- mean((t(vT[, selT]) - vtT[selT])^2) 
  # MSE Zero Coefficients Theta
  MsezT  <- mean((t(vT[, !selT]) - vtT[!selT])^2)
  # MSE Theta
  MseT <- mean((t(vT) - vtT)^2)
  
  # MSE Global and Local Non-Zero Coefficients B
  vB      <- B
  dim(vB) <- c(S, V * P)
  MseNzB  <- mean((t(vB[, c(tgB) == 1]) - tB[tgB == 1])^2)
  # MSE Local Zero
  vNzB      <- B[,,tg == 1]
  dim(vNzB) <- c(S, V * sum(tg == 1)) 
  MselzB    <- mean((vNzB[, tgB[, tg == 1] == 0])^2)
  # MSE Global Zero
  vNzB      <- B[,,tg == 0] 
  MsegzB    <- mean(vNzB^2)
  # MSE B
  MseB <- mean((t(vB) - c(tB))^2)
  
  # Global MSE
  Mse <- mean(c((t(vT) - vtT)^2, (t(vB) - c(tB))^2)) 
  
  # Coverage
  # Coverage Non-Zero Theta
  uppT   <- apply(vT, 2, quantile, p = 0.975)
  lowT   <- apply(vT, 2, quantile, p = 0.025)
  covNzT <- mean((lowT[selT] <= vtT[selT]) * (uppT[selT] >= vtT[selT]))
  # Coverage Zero Theta
  covzT <- mean((lowT[!selT] <= vtT[!selT]) * (uppT[!selT] >= vtT[!selT]))
  # Coverage Theta
  covT <- mean((lowT <= vtT) * (uppT >= vtT))
  
  # Coverage Non-Zero B
  uppB   <- apply(vB, 2, quantile, p = 0.975)
  lowB   <- apply(vB, 2, quantile, p = 0.025)
  covNzB <- mean((lowB[c(tgB == 1)] <= tB[tgB == 1]) * (uppB[c(tgB == 1)] >= tB[tgB == 1]))
  # Coverage Local Zero B
  vNzB      <- B[,,tg == 1]
  dim(vNzB) <- c(S, V * sum(tg == 1))
  upplzB   <- apply(vNzB, 2, quantile, p = 0.975)
  lowlzB   <- apply(vNzB, 2, quantile, p = 0.025)
  covlzB <- mean((lowlzB[c(tgB[, tg == 1] == 0)] <= tB[, tg == 1][tgB[, tg == 1] == 0]) * (upplzB[c(tgB[, tg == 1] == 0)] >= tB[, tg == 1][tgB[, tg == 1] == 0]))
  # Coverage Global Zero B
  vNzB      <- B[,,tg == 0]
  dim(vNzB) <- c(S, V * sum(tg == 0))
  uppgzB    <- apply(vNzB, 2, quantile, p = 0.975)
  lowgzB    <- apply(vNzB, 2, quantile, p = 0.025)
  covgzB    <- mean((lowgzB <= c(tB[, tg == 0]) * (uppgzB >= tB[, tg == 0])))
  # Coverage B
  covB    <- mean((lowB <= c(tB) * (uppB >= tB)))
  
  # Global Coverage
  cov <- mean(c((lowB <= c(tB) * (uppB >= tB)), (lowT <= vtT) * (uppT >= vtT)))
  
  # Length
  # Length Non-Zero Theta
  lenNzT <- mean(uppT[selT] - lowT[selT])
  # Coverage Zero Theta
  lenzT  <- mean(uppT[!selT] - lowT[!selT])
  # Coverage Theta
  lenT   <- mean(uppT - lowT)
  
  # Coverage Non-Zero B
  uppB   <- apply(vB, 2, quantile, p = 0.975)
  lowB   <- apply(vB, 2, quantile, p = 0.025)
  lenNzB <- mean(uppB[c(tgB == 1)] - lowB[c(tgB == 1)])
  # Coverage Local Zero B
  vNzB      <- B[,,tg == 1]
  dim(vNzB) <- c(S, V * sum(tg == 1))
  upplzB    <- apply(vNzB, 2, quantile, p = 0.975)
  lowlzB    <- apply(vNzB, 2, quantile, p = 0.025)
  lenlzB    <- mean(upplzB[c(tgB[, tg == 1] == 0)] - lowlzB[c(tgB[, tg == 1] == 0)])
  # Coverage Global Zero B
  vNzB      <- B[,,tg == 0]
  dim(vNzB) <- c(S, V * sum(tg == 0))
  uppgzB    <- apply(vNzB, 2, quantile, p = 0.975)
  lowgzB    <- apply(vNzB, 2, quantile, p = 0.025)
  lengzB    <- mean(uppgzB - lowgzB)
  # Coverage B
  lenB    <- mean(uppB - lowB)
  
  # Global Length
  len <- mean(c(uppB - lowB, uppT - lowT))
  
  
  # Return Stats
  return(list(TPR    = TPR,
              TNR    = TNR,
              MseNzT = MseNzT,
              MsezT  = MsezT,
              MseT   = MseT,
              MseNzB = MseNzB,
              MselzB = MselzB,
              MsegzB = MsegzB,
              MseB   = MseB,
              Mse    = Mse,
              covNzT = covNzT,
              covzT  = covzT,
              covT   = covT,
              covNzB = covNzB,
              covlzB = covlzB,
              covgzB = covgzB,
              covB   = covB,
              cov    = cov,
              lenNzT = lenNzT,
              lenzT  = lenzT,
              lenT   = lenT,
              lenNzB = lenNzB,
              lenlzB = lenlzB,
              lengzB = lengzB,
              lenB   = lenB,
              len    = len))
}