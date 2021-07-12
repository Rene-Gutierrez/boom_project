source('./data_generator4.R')
source('./fast_sampler.R')
source('./gboom_iterator.R')
source('./group_iterator.R')
source('./gboom_sampler.R')
source('./fast_horseshoe.R')
source('./horseshoe_sampler.R')
source('./gboom_stats2.R')
source('./pre_stats2.R')
source('./pre_data_gen.R')
source('./org_pre.R')
source('./org_sta.R')
source('./sta_lat.R')
source('./sim_boom.R')
source('./mcp_estimation.R')
source('./tra_sta.R')
source('./tra_pre_stats.R')
library(mclust)
library(ncvreg)

# For Replicability
set.seed(29062021)

# print(paste0("1st Case"))
# sim_boom(PP = 20,
#          VV = 20,
#          pB = 0.4,
#          pT = 0.25,
#          cB = 1,
#          cT = 3,
#          ss = 1,
#          nn = 150,
#          S  = 10000,
#          bI = 5001,
#          R  = 0,
#          N  = 100,
#          maxCoe = 100,
#          m      = 10)
# 
# print(paste0("2nd Case"))
# sim_boom(PP = 20,
#          VV = 15,
#          pB = 0.4,
#          pT = 0.25,
#          cB = 1,
#          cT = 3,
#          ss = 1,
#          nn = 150,
#          S  = 10000,
#          bI = 5001,
#          R  = 0,
#          N  = 100,
#          maxCoe = 100,
#          m      = 10)
# 
# print(paste0("3rd Case"))
# sim_boom(PP = 20,
#          VV = 10,
#          pB = 0.4,
#          pT = 0.25,
#          cB = 1,
#          cT = 3,
#          ss = 1,
#          nn = 150,
#          S  = 10000,
#          bI = 5001,
#          R  = 0,
#          N  = 100,
#          maxCoe = 100,
#          m      = 10)

print(paste0("4th Case"))
sim_boom(PP = 20,
         VV = 5,
         pB = 0.4,
         pT = 0.25,
         cB = 1,
         cT = 3,
         ss = 1,
         nn = 150,
         S  = 10000,
         bI = 5001,
         R  = 0,
         N  = 100,
         maxCoe = 100,
         m      = 10)