org_sta <- function(sta, filNam){
  # Number of Procedures
  M <- length(sta)
  N <- length(sta[[1]])
  
  # Statistics
  TR     <- matrix(data = NA, nrow = M, ncol = 2)
  TRSE   <- matrix(data = NA, nrow = M, ncol = 2)
  TRLab  <- c("TPR", "TNR")
  MSE    <- matrix(data = NA, nrow = M, ncol = 8)
  MSESE  <- matrix(data = NA, nrow = M, ncol = 8)
  MSELab <- c("Active", "Inactive", "All", "Active", "L. Inactive", "G. Inactive", "All", "All")
  COV    <- matrix(data = NA, nrow = M, ncol = 8)
  COVSE  <- matrix(data = NA, nrow = M, ncol = 8)
  COVLab <- c("Active", "Inactive", "All", "Active", "L. Inactive", "G. Inactive", "All", "All")
  LEN    <- matrix(data = NA, nrow = M, ncol = 8)
  LENSE  <- matrix(data = NA, nrow = M, ncol = 8)
  LENLab <- c("Active", "Inactive", "All", "Active", "L. Inactive", "G. Inactive", "All", "All")
  ind_sta <- c()
  ROC    <- array(data = NA, dim = c(M, N, 2, dim(sta[[1]][[1]]$ROC)[2]))
  for(m in 1:M){
    cur_sta <- sta[[m]]
    if(m != 3){
      # ROC
      for(n in 1:N){
        ROC[m, n, , ] <- cur_sta[[n]]$ROC
      }
      # Convergence Check
      mR <- numeric(N)
      for(n in 1:N){
        mR[n] <- cur_sta[[n]]$mR
      }
      # Effective Sample Size
      eff <- numeric(N)
      for(n in 1:N){
        eff[n] <- cur_sta[[n]]$avgEff
      }
      # Number of Coefficients
      nC <- numeric(N)
      for(n in 1:N){
        nC[n] <- cur_sta[[n]]$nC
      }
    }
    # TPR and FPR
    TPR     <- numeric(N)
    TNR     <- numeric(N)
    for(n in 1:N){
      TPR[n] <- cur_sta[[n]]$TPR
      TNR[n] <- cur_sta[[n]]$TNR
    }
    TR[m, 1]   <- mean(TPR)
    TR[m, 2]   <- mean(TNR)
    TRSE[m, 1] <- sd(TPR)
    TRSE[m, 2] <- sd(TNR)
    
    # MSE
    MseNzT <- numeric(N)
    MsezT  <- numeric(N)
    MseT   <- numeric(N)
    MseNzB <- numeric(N)
    MselzB <- numeric(N)
    MsegzB <- numeric(N)
    MseB   <- numeric(N)
    Mse    <- numeric(N)
    for(n in 1:N){
      MseNzT[n] <- cur_sta[[n]]$MseNzT
      MsezT[n]  <- cur_sta[[n]]$MsezT
      MseT[n]   <- cur_sta[[n]]$MseT
      MseNzB[n] <- cur_sta[[n]]$MseNzB
      MselzB[n] <- cur_sta[[n]]$MselzB
      MsegzB[n] <- cur_sta[[n]]$MsegzB
      MseB[n]   <- cur_sta[[n]]$MseB
      Mse[n]    <- cur_sta[[n]]$Mse
    }
    MSE[m, 1]   <- mean(MseNzT)
    MSE[m, 2]   <- mean(MsezT)
    MSE[m, 3]   <- mean(MseT)
    MSE[m, 4]   <- mean(MseNzB)
    MSE[m, 5]   <- mean(MselzB)
    MSE[m, 6]   <- mean(MsegzB)
    MSE[m, 7]   <- mean(MseB)
    MSE[m, 8]   <- mean(Mse)
    MSESE[m, 1] <- sd(MseNzT)
    MSESE[m, 2] <- sd(MsezT)
    MSESE[m, 3] <- sd(MseT)
    MSESE[m, 4] <- sd(MseNzB)
    MSESE[m, 5] <- sd(MselzB)
    MSESE[m, 6] <- sd(MsegzB)
    MSESE[m, 7] <- sd(MseB)
    MSESE[m, 8] <- sd(Mse)
    
    if(m != 3){
      # Coverage
      covNzT <- numeric(N)
      covzT  <- numeric(N)
      covT   <- numeric(N)
      covNzB <- numeric(N)
      covlzB <- numeric(N)
      covgzB <- numeric(N)
      covB   <- numeric(N)
      cov    <- numeric(N)
      for(n in 1:N){
        covNzT[n] <- cur_sta[[n]]$covNzT
        covzT[n]  <- cur_sta[[n]]$covzT
        covT[n]   <- cur_sta[[n]]$covT
        covNzB[n] <- cur_sta[[n]]$covNzB
        covlzB[n] <- cur_sta[[n]]$covlzB
        covgzB[n] <- cur_sta[[n]]$covgzB
        covB[n]   <- cur_sta[[n]]$covB
        cov[n]    <- cur_sta[[n]]$cov
      }
      COV[m, 1]   <- mean(covNzT)
      COV[m, 2]   <- mean(covzT)
      COV[m, 3]   <- mean(covT)
      COV[m, 4]   <- mean(covNzB)
      COV[m, 5]   <- mean(covlzB)
      COV[m, 6]   <- mean(covgzB)
      COV[m, 7]   <- mean(covB)
      COV[m, 8]   <- mean(cov)
      COVSE[m, 1] <- sd(covNzT)
      COVSE[m, 2] <- sd(covzT)
      COVSE[m, 3] <- sd(covT)
      COVSE[m, 4] <- sd(covNzB)
      COVSE[m, 5] <- sd(covlzB)
      COVSE[m, 6] <- sd(covgzB)
      COVSE[m, 7] <- sd(covB)
      COVSE[m, 8] <- sd(cov)
      
      # Length
      lenNzT <- numeric(N)
      lenzT  <- numeric(N)
      lenT   <- numeric(N)
      lenNzB <- numeric(N)
      lenlzB <- numeric(N)
      lengzB <- numeric(N)
      lenB   <- numeric(N)
      len    <- numeric(N)
      for(n in 1:N){
        lenNzT[n] <- cur_sta[[n]]$lenNzT
        lenzT[n]  <- cur_sta[[n]]$lenzT
        lenT[n]   <- cur_sta[[n]]$lenT
        lenNzB[n] <- cur_sta[[n]]$lenNzB
        lenlzB[n] <- cur_sta[[n]]$lenlzB
        lengzB[n] <- cur_sta[[n]]$lengzB
        lenB[n]   <- cur_sta[[n]]$lenB
        len[n]    <- cur_sta[[n]]$len
      }
      LEN[m, 1]   <- mean(lenNzT)
      LEN[m, 2]   <- mean(lenzT)
      LEN[m, 3]   <- mean(lenT)
      LEN[m, 4]   <- mean(lenNzB)
      LEN[m, 5]   <- mean(lenlzB)
      LEN[m, 6]   <- mean(lengzB)
      LEN[m, 7]   <- mean(lenB)
      LEN[m, 8]   <- mean(len)
      LENSE[m, 1] <- sd(lenNzT)
      LENSE[m, 2] <- sd(lenzT)
      LENSE[m, 3] <- sd(lenT)
      LENSE[m, 4] <- sd(lenNzB)
      LENSE[m, 5] <- sd(lenlzB)
      LENSE[m, 6] <- sd(lengzB)
      LENSE[m, 7] <- sd(lenB)
      LENSE[m, 8] <- sd(len)
    }
    ind_sta <- cbind(ind_sta,
                     TPR, TNR,
                     Mse, MseNzT, MsezT, MseT, MseNzB, MselzB, MsegzB, MseB,
                     cov, covNzT, covzT, covT, covNzB, covlzB, covgzB, covB,
                     len, lenNzT, lenzT, lenT, lenNzB, lenlzB, lengzB, lenB,
                     mR, eff)
  }
  ind_sta <- cbind(ind_sta, nC)
  write.table(ind_sta, file = paste0(filNam, '.txt'), quote = FALSE)
  saveRDS(ind_sta, file = paste0(filNam, '.rds'))
  saveRDS(ROC, file = paste0("ROC_",filNam, '.rds'))
  # Returns
  return(list(TR     = TR,
              TRSE   = TRSE,
              TRLab  = TRLab,
              MSE    = MSE,
              MSESE  = MSESE,
              MSELab = MSELab,
              COV    = COV,
              COVSE  = COVSE,
              COVLab = COVLab,
              LEN    = LEN,
              LENSE  = LENSE,
              LENLab = LENLab))
}