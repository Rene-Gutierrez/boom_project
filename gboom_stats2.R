gboom_stats2 <- function(tTheta,
                        tB,
                        tg,
                        tgB,
                        Theta,
                        B,
                        g       = NA,
                        flag1   = TRUE,
                        m,
                        lag.max = 50,
                        nC,
                        cO = 1 - seq(0,1, length.out = 101)){
  # Dimensions
  S <- dim(B)[1]
  V <- dim(B)[2]
  P <- dim(B)[3]
  
  if(flag1){
    # True Positive Rate
    mg  <- apply(X = g, MARGIN = 2, FUN = median)
    TPR <- mean(mg[tg == 1])
    TNR <- 1 - mean(mg[tg == 0])
    # ROC
    ncO <- length(cO)
    ROC <- matrix(data = NA, nrow = 2, ncol = ncO)
    for(i in 1:ncO){
      ROC[1, i] <- mean(colMeans(g[, tg == 1]) >  cO[i])
      ROC[2, i] <- mean(colMeans(g[, tg == 0]) >= cO[i])
    }
  } else {
    # Cluster T
    cluRes <- mclust::Mclust(data = abs(c(apply(Theta, c(2,3), mean))), G = 1:2, verbose = FALSE)
    selT   <- colSums(matrix(cluRes$classification, P, P) == which.max(cluRes$parameters$mean)) > 0
    # Cluster B
    cluRes <- mclust::Mclust(data = abs(c(apply(B, c(2,3), mean))), G = 1:2, verbose = FALSE)
    selB   <- colSums(matrix(cluRes$classification, V, P) == which.max(cluRes$parameters$mean)) > 0
    # Final Selection
    g      <- ((selT + selB) > 0) + 0
    TPR <- mean(g[tg == 1])
    TNR <- 1 - mean(g[tg == 0])
    # ROC
    ncO <- length(cO)
    ROC <- matrix(data = NA, nrow = 2, ncol = ncO)
  }
  
  
  # MSE Non-Zero Coefficients Theta
  vT      <- Theta
  dim(vT) <- c(S, P * P)
  vT      <- vT[, lower.tri(Theta[1,,])]
  vT      <- apply(X = vT, MARGIN = 2, FUN = median)
  vtT     <- tTheta[lower.tri(tTheta)]
  selT    <- ((tg %*% t(tg)) == 1)
  selT    <- selT[lower.tri(selT)]
  MseNzT  <- mean((t(vT[selT]) - vtT[selT])^2) 
  # MSE Zero Coefficients Theta
  MsezT  <- mean((t(vT[!selT]) - vtT[!selT])^2)
  # MSE Theta
  MseT <- mean((t(vT) - vtT)^2)
  
  # MSE Global and Local Non-Zero Coefficients B
  vB      <- B
  dim(vB) <- c(S, V * P)
  vB      <- apply(X = vB, MARGIN = 2, FUN = median)
  MseNzB  <- mean((t(vB[c(tgB) == 1]) - tB[tgB == 1])^2)
  # MSE Local Zero
  vNzB      <- B[,,tg == 1]
  dim(vNzB) <- c(S, V * sum(tg == 1))
  vNzB      <- apply(X = vNzB, MARGIN = 2, FUN = median)
  MselzB    <- mean((vNzB[tgB[, tg == 1] == 0])^2)
  # MSE Global Zero
  vNzB      <- B[,,tg == 0]
  vNZB      <- apply(X = vNzB, MARGIN = c(2, 3), FUN = median)
  MsegzB    <- mean(vNzB^2)
  # MSE B
  MseB <- mean((vB - c(tB))^2)
  
  # Global MSE
  Mse <- mean(c((vT - vtT)^2, (vB - c(tB))^2)) 
  
  # Coverage
  vT      <- Theta
  dim(vT) <- c(S, P * P)
  vT      <- vT[, lower.tri(Theta[1,,])]
  # Coverage Non-Zero Theta
  uppT   <- apply(vT, 2, quantile, p = 0.975)
  lowT   <- apply(vT, 2, quantile, p = 0.025)
  covNzT <- mean((lowT[selT] <= vtT[selT]) * (uppT[selT] >= vtT[selT]))
  # Coverage Zero Theta
  covzT <- mean((lowT[!selT] <= vtT[!selT]) * (uppT[!selT] >= vtT[!selT]))
  # Coverage Theta
  covT <- mean((lowT <= vtT) * (uppT >= vtT))
  
  # Coverage Non-Zero B
  vB      <- B
  dim(vB) <- c(S, V * P)
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
  
  # Convergence Evaluation B
  n   <- S / m
  ind <- c(t(matrix(rep(1:m, n), nrow = m, ncol = n)))
  mB  <- matrix(data = NA, nrow = m, ncol = P * V)
  sB  <- matrix(data = NA, nrow = m, ncol = P * V)
  for(i in 1:m){
    miniB   <- vB[ind == i,]
    mB[i, ] <- colMeans(miniB)
    sB[i, ] <- apply(miniB, 2, var)
  }
  B_B      <- n * apply(mB, 2, var)
  W_B      <- colMeans(sB)
  varHat_B <- (n - 1) * W_B / n + 1 * B_B / n
  R_B      <- sqrt(varHat_B / W_B)
  mR_B     <- sqrt(((n - 1) * sum(W_B) / n + 1 * sum(B_B) / n) / sum(W_B))
  
  # Convergence Evaluation Theta
  mT  <- matrix(data = NA, nrow = m, ncol = P * (P -1) / 2)
  sT  <- matrix(data = NA, nrow = m, ncol = P * (P -1) / 2)
  for(i in 1:m){
    miniT   <- vT[ind == i,]
    mT[i, ] <- colMeans(miniT)
    sT[i, ] <- apply(miniT, 2, var)
  }
  B_T      <- n * apply(mT, 2, var)
  W_T      <- colMeans(sT)
  varHat_T <- (n - 1) * W_T / n + 1 * B_T / n
  R_T      <- sqrt(varHat_T / W_T)
  mR_T     <- sqrt(((n - 1) * sum(W_T) / n + 1 * sum(B_T) / n) / sum(W_T))
  
  mR <- (mR_B * (P * V) + mR_T * P * (P - 1) / 2) / (P * V + P * (P - 1) / 2)
  
  # Effective Number of Samples
  vT      <- Theta
  dim(vT) <- c(S, P * P)
  vT      <- vT[, lower.tri(Theta[1,,])]
  selT    <- ((tg %*% t(tg)) == 1)
  selT    <- selT[lower.tri(selT)]
  selT    <- vT[, selT]
  vB      <- B
  dim(vB) <- c(S, V * P)
  selB    <- vB[, c(tgB == 1)]
  selBT   <- cbind(selT, selB)
  ncolS   <- dim(selBT)[2]
  eff     <- numeric(length = ncolS)
  for(i in 1:ncolS){
    acfS   <- c(acf(x = selBT[, i], lag.max = lag.max, plot = FALSE)$acf)
    cutOff <- which.max((acfS[1:lag.max] < 0) * (acfS[2:(lag.max + 1)] < 0))
    if(length(cutOff) == 0){
      cutOff <- 1
    }
    sumACF <- sum(acfS[2:(cutOff + 1)])
    effNum <- 1 + 2 * (sumACF)
    eff[i] <- 1 / effNum 
  }
  
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
              len    = len,
              varHat_B = varHat_B,
              W_B    = W_B,
              B_B    = B_B,
              mR_B   = mR_B,
              varHat_T = varHat_T,
              W_T    = W_T,
              B_T    = B_T,
              mR_T   = mR_T,
              mR     = mR,
              eff    = eff,
              avgEff = mean(eff),
              nC     = nC,
              ROC    = ROC))
}