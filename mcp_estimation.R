mcp_estimation <- function(y, G, A){
  # Problem Dimensions
  n <- dim(G)[1]  # Number of Observations
  V <- dim(G)[2]  # Voxel Size
  P <- dim(G)[3]  # Number of ROI's
  
  # Auxiliary Variables
  GX                <- G
  dim(GX)           <- c(n, V * P)
  AX                <- A
  dim(AX)           <- c(n, P * P)
  AX                <- AX[, lower.tri(A[1,,])]
  cvI               <- matrix(data = sample(x = 1:n, size = n), nrow = n / 10, ncol = 10)
  err               <- matrix(data = NA, nrow = 10, ncol = 100)
  
  for(i in 1:10){
    val <- cvI[, i]
    out <- ncvreg::ncvreg(X = cbind(GX[-val,], AX[-val,]),
                          y = y[-val])$beta[-1, ]
    OSe <- colMeans((cbind(GX[val,], AX[val,]) %*% out - y[val])^2)
    err[i, ] <- OSe
  }
  err <- colMeans(err)
  lam <- which.min(err)
  
  coe <- ncvreg::ncvreg(X = cbind(GX, AX),
                        y = y)$beta[-1, lam]
  
  return(list(GX  = GX,
              AX  = AX,
              y   = y,
              cvI = cvI,
              err = err,
              lam = lam,
              coe = coe))
}