source('./data_generator.R')
source('./fast_sampler.R')
source('./fast_horseshoe.R')
source('./boom_sampler.R')

data_out <- data_generator(P = 20,   # Number of Regions of Interest
                           V = 25,   # Number of Voxels per region
                           R = 2,    # Depth of the Connectivity Matrix
                           a = 0,    # Signal Connectivity Matrix
                           W = 2,    # Active Voxels per Region
                           c = 10,   # Signal per Voxel
                           Q = 10,   # Active Regions of Interest
                           s = 1,    # Regression Noise
                           n = 100)  # Number of Samples

X <- data_out$G
dim(X) <- c(100, 25 * 20)
y <- data_out$y

horse_out <- horseshoe::horseshoe(y = y, X = X, method.tau = "halfCauchy", method.sigma = "Jeffreys", burn = 0, nmc = 10000)

out <- boom_sampler_old(A = data_out$A,
                    G = data_out$G,
                    y = data_out$y,
                    S = 10000)

burnIn <- 5000
mcmc   <- 10000
ymin   <- min(out$B[(burnIn + 1):mcmc,], horse_out$BetaSamples[, (burnIn + 1):mcmc])
ymax   <- max(out$B[(burnIn + 1):mcmc,], horse_out$BetaSamples[, (burnIn + 1):mcmc])
par(mfcol = c(3, 2))
plot(out$B[(burnIn + 1):mcmc, 101], type = 'l', ylim = c(ymin, ymax))
abline(h = 1, col = "red", lwd = 2)
plot(out$B[(burnIn + 1):mcmc, 102], type = 'l', ylim = c(ymin, ymax))
abline(h = 1, col = "red", lwd = 2)
plot(out$B[(burnIn + 1):mcmc, 103], type = 'l', ylim = c(ymin, ymax))
abline(h = 0, col = "red", lwd = 2)
plot(horse_out$BetaSamples[101, (burnIn + 1):mcmc], type = 'l', ylim = c(ymin, ymax))
abline(h = 1, col = "red", lwd = 2)
plot(horse_out$BetaSamples[102, (burnIn + 1):mcmc], type = 'l', ylim = c(ymin, ymax))
abline(h = 1, col = "red", lwd = 2)
plot(horse_out$BetaSamples[103, (burnIn + 1):mcmc], type = 'l', ylim = c(ymin, ymax))
abline(h = 0, col = "red", lwd = 2)
plot(log(out$l2[(burnIn + 1):mcmc, 1]), type = 'l')
plot(log(out$l2[(burnIn + 1):mcmc, 2]), type = 'l')
plot(log(out$l2[(burnIn + 1):mcmc, 3]), type = 'l')
plot(log(out$v[(burnIn + 1):mcmc, 1]), type = 'l')
plot(log(out$v[(burnIn + 1):mcmc, 2]), type = 'l')
plot(log(out$v[(burnIn + 1):mcmc, 3]), type = 'l')
par(mfcol = c(2, 2))
plot(log(out$s2[(burnIn + 1):mcmc]), type = 'l')
plot(log(out$t2[(burnIn + 1):mcmc]), type = 'l')
plot(log(horse_out$Sigma2Samples[(burnIn + 1):mcmc]), type = 'l')
plot(log(horse_out$TauSamples[(burnIn + 1):mcmc]), type = 'l')

