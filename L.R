### Prior Correlation Matrix

L <- function(M, g, a, m2, M2){
  Bpart     <- c(t((1 - g) * a + g * t(M)))
  Thetapart <- g %*% t(g)
  Thetapart <- Thetapart[lower.tri(Thetapart)]
  Thetapart <- (1 - Thetapart) * m2 + Thetapart * M2
  return(c(Bpart, Thetapart))
}