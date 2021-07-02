org_pre <- function(pre, filNam){
  # Number of Procedures
  M <- length(pre)
  N <- length(pre[[1]])
  
  # Statistics
  PRE    <- matrix(data = NA, nrow = M, ncol = 4)
  PRESE  <- matrix(data = NA, nrow = M, ncol = 4)
  PRELab <- c("MSE", "MSE \\%", "Cov", "Len")
  
  ind_sta <- c()
  for(m in 1:M){
    cur_sta <- pre[[m]]
    # Stats
    MSE <- numeric(N)
    PMS <- numeric(N)
    COV <- numeric(N)
    LEN <- numeric(N)
      for(n in 1:N){
      MSE[n] <- cur_sta[[n]]$mse
      PMS[n] <- cur_sta[[n]]$pmse
      COV[n] <- cur_sta[[n]]$cov
      LEN[n] <- cur_sta[[n]]$len
      }
    PRE[m, 1] <- mean(MSE)
    PRE[m, 2] <- mean(PMS) * 100
    PRE[m, 3] <- mean(COV)
    PRE[m, 4] <- mean(LEN)
    PRESE[m, 1] <- sd(MSE)
    PRESE[m, 2] <- sd(PMS) * 100
    PRESE[m, 3] <- sd(COV)
    PRESE[m, 4] <- sd(LEN)
    ind_sta <- cbind(ind_sta, MSE, PMS, COV, LEN)
  }
  write.table(ind_sta, file = filNam, quote = FALSE)
  # Returns the Stats
  return(list(PRE   = PRE,
              PRESE = PRESE,
              PRELab = PRELab))
}