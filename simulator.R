source('./data_generator2.R')
source('./fast_sampler.R')
source('./gboom_iterator.R')
source('./group_iterator.R')
source('./gboom_sampler.R')
source('./fast_horseshoe.R')
source('./horseshoe_sampler.R')
source('./gboom_stats.R')
source('./pre_stats.R')
source('./pre_data_gen.R')
source('./org_pre.R')
source('./org_sta.R')
source('./sta_lat.R')
source('./sim_boom.R')

# For Replicability
set.seed(29062021)

print(paste0("1st Case"))
sim_boom(PP = 20,
         VV = 15,
         pB = 0.2,
         pT = 0.2,
         cB = 3,
         cT = 1,
         ss = 1,
         nn = 100,
         S  = 2000,
         bI = 1001,
         R  = 0,
         N  = 20)

print(paste0("2nd Case"))
sim_boom(PP = 20,
         VV = 10,
         pB = 0.2,
         pT = 0.2,
         cB = 3,
         cT = 1,
         ss = 1,
         nn = 100,
         S  = 2000,
         bI = 1001,
         R  = 0,
         N  = 20)

print(paste0("3rd Case"))
sim_boom(PP = 20,
         VV = 5,
         pB = 0.2,
         pT = 0.2,
         cB = 3,
         cT = 1,
         ss = 1,
         nn = 100,
         S  = 2000,
         bI = 1001,
         R  = 0,
         N  = 20)

print(paste0("4th Case"))
sim_boom(PP = 20,
         VV = 15,
         pB = 0.2,
         pT = 0.3,
         cB = 3,
         cT = 1,
         ss = 1,
         nn = 100,
         S  = 2000,
         bI = 1001,
         R  = 0,
         N  = 20)

print(paste0("5th Case"))
sim_boom(PP = 20,
         VV = 10,
         pB = 0.2,
         pT = 0.3,
         cB = 3,
         cT = 1,
         ss = 1,
         nn = 100,
         S  = 2000,
         bI = 1001,
         R  = 0,
         N  = 20)

print(paste0("6th Case"))
sim_boom(PP = 20,
         VV = 5,
         pB = 0.2,
         pT = 0.3,
         cB = 3,
         cT = 1,
         ss = 1,
         nn = 100,
         S  = 2000,
         bI = 1001,
         R  = 0,
         N  = 20)