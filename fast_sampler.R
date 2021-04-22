### Fast Sampler

fast_sampler <- function(Phi, D, a){
  # Dimensions
  n <- dim(Phi)[1]
  p <- dim(Phi)[2]
  
  # Step 1
  u <- rnorm(n = p, mean = 0, sd = sqrt(D))
  d <- rnorm(n = n, mean = 0, sd = 1)
  
  # Step 2
  v <- Phi %*% u + d
  
  # Step 3
  w <- solve(Phi %*%  (D * t(Phi)) + diag(n), a - v)
  
  # Step 4
  t <- u + D * t(Phi) %*% w
  
  # Returns the Sample
  return(t)
}