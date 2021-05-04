source('./data_generator.R')
source('./fast_sampler.R')
source('./fast_spike_horseshoe.R')
source('./boom_sampler.R')

data_out <- data_generator(P = 20,   # Number of Regions of Interest
                           V = 10,   # Number of Voxels per region
                           R = 2,    # Depth of the Connectivity Matrix
                           a = 3,   # Signal Connectivity Matrix
                           W = 2,    # Active Voxels per Region
                           c = 1,    # Signal per Voxel
                           Q = 5,    # Active Regions of Interest
                           s = 1,    # Regression Noise
                           n = 100)  # Number of Samples

S <- 10000

out <- boom_sampler(A = data_out$A,
                    G = data_out$G,
                    y = data_out$y,
                    S = S)

sta <- 1
fin <- S

par(mfcol = c(2, 2))
# B
plot(out$B[sta:fin, 1, 1], type = 'l')
abline(h = 1, col = "red", lwd = 2)
plot(out$B[sta:fin, 2, 1], type = 'l')
abline(h = 1, col = "red", lwd = 2)
plot(out$B[sta:fin, 1, 6], type = 'l')
abline(h = 0, col = "red", lwd = 2)
plot(out$B[sta:fin, 1, 7], type = 'l')
abline(h = 0, col = "red", lwd = 2)
# l2
plot(log(out$l2[sta:fin, 1, 1]), type = 'l')
plot(log(out$l2[sta:fin, 2, 1]), type = 'l')
plot(log(out$l2[sta:fin, 1, 6]), type = 'l')
plot(log(out$l2[sta:fin, 1, 7]), type = 'l')
# v
plot(log(out$v[sta:fin, 1, 1]), type = 'l')
plot(log(out$v[sta:fin, 2, 1]), type = 'l')
plot(log(out$v[sta:fin, 1, 6]), type = 'l')
plot(log(out$v[sta:fin, 1, 7]), type = 'l')
# t2
plot(log(out$t2[sta:fin, 1]), type = 'l')
plot(log(out$t2[sta:fin, 2]), type = 'l')
plot(log(out$t2[sta:fin, 6]), type = 'l')
plot(log(out$t2[sta:fin, 7]), type = 'l')
# t2 * l2
plot(log(out$l2_t2[sta:fin, 1, 1]), type = 'l')
abline(h = log(0.01), col = "red", lwd = 2)
plot(log(out$l2_t2[sta:fin, 1, 2]), type = 'l')
abline(h = log(0.01), col = "red", lwd = 2)
plot(log(out$l2_t2[sta:fin, 1, 6]), type = 'l')
abline(h = log(0.01), col = "red", lwd = 2)
plot(log(out$l2_t2[sta:fin, 1, 7]), type = 'l')
abline(h = log(0.01), col = "red", lwd = 2)
# e2
plot(out$e2[sta:fin, 1], type = 'l')
plot(out$e2[sta:fin, 2], type = 'l')
plot(out$e2[sta:fin, 6], type = 'l')
plot(out$e2[sta:fin, 7], type = 'l')
# xi
plot(log(out$xi[sta:fin, 1]), type = 'l')
plot(log(out$xi[sta:fin, 2]), type = 'l')
plot(log(out$xi[sta:fin, 6]), type = 'l')
plot(log(out$xi[sta:fin, 7]), type = 'l')
# g
plot(out$g[sta:fin, 1], type = 'l')
plot(out$g[sta:fin, 2], type = 'l')
plot(out$g[sta:fin, 6], type = 'l')
plot(out$g[sta:fin, 7], type = 'l')
# pg
plot(out$pg[sta:fin, 1], type = 'l')
plot(out$pg[sta:fin, 2], type = 'l')
plot(out$pg[sta:fin, 6], type = 'l')
plot(out$pg[sta:fin, 7], type = 'l')
# s2
par(mfcol = c(1, 1))
plot(out$s2[sta:fin], type = 'l')
