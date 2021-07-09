sim_boom <- function(PP = 20,
                     VV = 10,
                     pB = 0.1,
                     pT = 0.2,
                     cB = 3,
                     cT = 1,
                     ss = 1,
                     nn = 100,
                     S  = 4000,
                     bI = 2001,
                     R  = 0,
                     N  = 100,
                     maxCoe = 100,
                     m      = 10){
  sta1 <- list()
  sta2 <- list()
  sta3 <- list()
  tim1 <- numeric(length = N)
  tim2 <- numeric(length = N)
  tim3 <- numeric(length = N)
  pre1 <- list()
  pre2 <- list()
  pre3 <- list()
  sele <- matrix(data = NA, nrow = N, ncol = S)
  unse <- matrix(data = NA, nrow = N, ncol = S)
  for(i in 1:N){
    # Tracker
    print(i)
    nC <- maxCoe + 1
    # Data Generation
    while (nC > maxCoe) {
      print(paste0("Try "))
      data_out <- data_generator4(P  = PP,
                                  V  = VV,
                                  pB = pB,
                                  pT = pT,
                                  cB = cB,
                                  cT = cT,
                                  s2 = ss,
                                  n  = nn)
      nC <- data_out$nC
      print(nC)
    }
    
    # Estimation
    print("G-Boom")
    tim1[i] <- Sys.time()
    out1 <- gboom_sampler(y  = data_out$y,
                          G  = data_out$G,
                          A  = data_out$A,
                          g  = rep(0, PP),
                          s2 = 1,
                          r  = 1 / 2,
                          S  = S,
                          R  = R)
    tim1[i]  <- Sys.time() - tim1[i]
    sele[i,] <- rowMeans(out1$g[, data_out$gT == 1])
    unse[i,] <- rowMeans(out1$g[, data_out$gT == 0])
    
    print("")
    print("Horseshoe")
    tim2[i] <- Sys.time()
    out2 <- horseshoe_sampler(A = data_out$A,
                              G = data_out$G,
                              y = data_out$y,
                              S = S)
    tim2[i]  <- Sys.time() - tim2[i]
    
    print("")
    print("MCP")
    tim3[i] <- Sys.time()
    out3 <- mcp_estimation(y = data_out$y,
                           G = data_out$G,
                           A = data_out$A)
    tim3[i]  <- Sys.time() - tim3[i]
    
    sta1[[i]] <- gboom_stats2(tTheta = data_out$Theta,
                             tB     = data_out$B,
                             tg     = data_out$gT,
                             tgB    = data_out$gB,
                             Theta  = out1$Theta[bI:S,,],
                             B      = out1$B[bI:S,,],
                             g      = out1$g[bI:S,],
                             m      = m,
                             nC     = data_out$nC)
    
    sta2[[i]] <- gboom_stats2(tTheta = data_out$Theta,
                             tB     = data_out$B,
                             tg     = data_out$gT,
                             tgB    = data_out$gB,
                             Theta  = out2$Theta[bI:S,,],
                             B      = out2$B[bI:S,,],
                             flag   = FALSE,
                             m      = m,
                             nC     = data_out$nC)
    
    sta3[[i]] <- tra_sta(tTheta = data_out$Theta,
                         tB     = data_out$B,
                         tg     = data_out$gT,
                         tgB    = data_out$gB,
                         coe    = out3$coe)
    
    # Predictive Stats
    # Predictions
    print("")
    print("Prediction")
    pre_dat <- pre_data_gen(B     = data_out$B,
                            Theta = data_out$Theta,
                            s2    = data_out$s2,
                            n     = 100)
    
    pre1[[i]] <- pre_stats2(ty    = pre_dat$y,
                           A     = pre_dat$A,
                           G     = pre_dat$G,
                           B     = out1$B[bI:S,,],
                           Theta = out1$Theta[bI:S,,])
    
    pre2[[i]] <- pre_stats2(ty    = pre_dat$y,
                           A     = pre_dat$A,
                           G     = pre_dat$G,
                           B     = out2$B[bI:S,,],
                           Theta = out2$Theta[bI:S,,])
    
    pre3[[i]] <- tra_pre_stats(ty = pre_dat$y,
                               A = pre_dat$A,
                               G = pre_dat$G,
                               coe = out3$coe)
  }
  
  suf <- paste0("_", pT, "_", pB, "_", VV, "_", PP, "_", nn, "_", cB, "_", cT, "_")
  saveRDS(cbind(tim1, tim2, tim3), file = paste0("tim_", suf, '.rds'))
  org <- org_sta(sta = list(sta1, sta2, sta3), filNam = paste0("per", suf))
  sta_lat(m = org$TR, se = org$TRSE, lab = org$TRLab, minmax = rep("max", 2), met = c("Boom", "Horseshoe", "MCP"), dig = 2, filNam = paste0("tr", suf))
  sta_lat(m = org$LEN, se = org$LENSE, lab = org$LENLab, minmax = rep("min", 8), met = c("Boom", "Horseshoe", "MCP"), dig = 3, filNam = paste0("len", suf))
  sta_lat(m = org$COV, se = org$COVSE, lab = org$MSELab, minmax = rep("max", 8), met = c("Boom", "Horseshoe", "MCP"), dig = 3, filNam = paste0("cov", suf))
  sta_lat(m = org$MSE, se = org$MSESE, lab = org$MSELab, minmax = rep("min", 8), met = c("Boom", "Horseshoe", "MCP"), dig = 2, filNam = paste0("mse", suf))
  
  pre <- org_pre(pre = list(pre1, pre2, pre3), filNam = paste0("per_pre", suf))
  sta_lat(m = pre$PRE, se = pre$PRESE, lab = pre$PRELab, minmax = c("min", "min", "max", "min"), met = c("Boom", "Horseshoe", "MCP"), dig = 4, filNam = paste0("pre", suf))
}