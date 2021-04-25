# 

# Inputs
P
n
R   <- 8
S   <- 100
tau <- 0.1

# Initializes u
u <- matrix(data = rnorm(R * P, mean = 0, sd = 0.1), nrow = R, ncol = P)

# Initializes lambda
lambda    <- rep(1, R)

# Error Tracking
e.y        <- numeric(length = S + 1)
e.Theta    <- numeric(length = S + 1)
e.Theta[1] <- sum(abs(r.Theta - t(u) %*% (lambda * u)))
temp       <- A
dim(temp)  <- c(n, P * P)
e.y[1]     <- mean(abs(y - rowSums(t(t(temp) * c(t(u) %*% (lambda * u)))) / 2))
start_time <- Sys.time()
for(s in 1:S){
  #  Estimates each u_p
  for(p in 1:P){
    ## Creates the Auxiliary Response Variable
    temp1      <- A[, -p, -p]
    dim(temp1) <- c(n, (P - 1) * (P - 1))
    temp2      <- t(u[, -p]) %*% (lambda * u[, -p])
    dim(temp2) <- (P - 1) * (P - 1)
    a.y        <- y - temp1 %*% temp2 / 2

    ## Creates the Auxiliary Covariate
    a.X <- A[, p, -p] %*% t(lambda * u[,-p])
    
    ## Point Estimates
    u[, p] <- solve( t(a.X) %*% a.X, t(a.X) %*% a.y)
    
    ## Error Tracking
    e.Theta[s + 1] <- sum(abs(r.Theta - t(u) %*% (lambda * u)))
    temp           <- A
    dim(temp)      <- c(n, P * P)
    e.y[s + 1]     <- mean(abs(y - rowSums(t(t(temp) * c(t(u) %*% (lambda * u)))) / 2))
  }
  
  # Estimates Lambda
  ## Creates an Auxiliary Covariate
  temp1      <- A
  dim(temp1) <- c(n, P * P)
  for(r in 1:R){
    temp2      <- u[r,] %*% t(u[r,])
    dim(temp2) <- P * P
    a.X        <- temp1 %*% temp2 / 2
    
    ## Creates an Auxiliary Response Variable
    temp3      <- matrix(data = u[-r,], nrow = R - 1, ncol = P)
    temp2      <- t(temp3) %*% temp3
    dim(temp2) <- P * P
    a.y        <- y - temp1 %*% temp2 / 2
    
    ## Point Estimates
    prob1 <- sum(dnorm(x = a.y, mean =  a.X, sd = tau, log = TRUE))
    prob2 <- sum(dnorm(x = a.y, mean =    0, sd = tau, log = TRUE))
    prob3 <- sum(dnorm(x = a.y, mean = -a.X, sd = tau, log = TRUE))
    pro1  <- prob1 - max(prob1, prob2, prob3)
    pro2  <- prob2 - max(prob1, prob2, prob3)
    pro3  <- prob3 - max(prob1, prob2, prob3)
    pro1  <- exp(pro1)
    pro2  <- exp(pro2)
    pro3  <- exp(pro3)
    pr1   <- pro1 / (pro1 + pro2 + pro3)
    pr2   <- pro2 / (pro1 + pro2 + pro3)
    pr3   <- pro3 / (pro1 + pro2 + pro3)
    lambda[r] <- c(1, 0, -1)[which.max(c(pr1, pr2, pr3))]
  }
}

# Plots Errors
plot(e.Theta, type = 'l', ylim = c(0, max(20)))
plot(e.y, type = 'l', ylim = c(0, max(5)))
end_time <- Sys.time()
print(end_time - start_time)