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
set.seed(20062021)

print(paste0("1st Case"))
sim_boom(pT = 0.2,
         VV = 5)

print(paste0("2nd Case"))
sim_boom(pT = 0.3,
         VV = 5)

print(paste0("3rd Case"))
sim_boom(pT = 0.2,
         VV = 10)

print(paste0("4th Case"))
sim_boom(pT = 0.3,
         VV = 10)

print(paste0("5th Case"))
sim_boom(pT = 0.2,
         VV = 15)

print(paste0("6th Case"))
sim_boom(pT = 0.3,
         VV = 15)