##### type I error and power analysis 

# data, functions and libraries -------------------------------------------
library(parallel)
source(here::here("R", "functions", "autoregressive_GN_01-09.R"))

# rho 0.0001 envir FALSE --------------------------------------------------


runs<- 1:1000 # number of tests
Nspp<- 50 
Ncomm<- 50
power<- 0.0001 # grafen parameter
perm<- 1000 # number of permutations in taxa shuffle null model
envir<- FALSE # environment
u<- 5
binary<- TRUE
test<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
}) # simulating trees - equal to the number of tests

ncor<- detectCores() # setting number of cores
CL<- makeCluster(ncor, setup_timeout = 0.5) 
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

system.time(res_rho.0001<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)) # run analysis
stopCluster(CL)
saveRDS(res_rho.0001, file = here::here("output", "Result_rho0001.rds"))



# rho - 0.8 envir FALSE ----------------------------------------------------

runs<- 1:1000 # number of tests
Nspp<- 50 
Ncomm<- 50
power<- 0.8 # grafen parameter
perm<- 1000 # number of permutations in taxa shuffle null model
envir<- FALSE # environment
u<- 5
binary<- TRUE
test<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
}) # simulating trees - equal to the number of tests

ncor<- detectCores() # setting number of cores
CL<- makeCluster(ncor, setup_timeout = 0.5) 
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

system.time(res_rho.8<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
  )) # run analysis
stopCluster(CL)
saveRDS(res_rho.8, file = here::here("output", "Result_rho08.rds"))

# rho 1.2 envir FALSE -------------------------------------------------

library(parallel)
runs<- 1:1000
Nspp<- 50
Ncomm<- 50
power<- 1.2
perm<- 1000
envir<- FALSE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
})
ncor<- detectCores()
CL<- makeCluster(ncor, setup_timeout = 0.5)
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho1.2<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)
saveRDS(res_rho1.2, file = here::here("output", "Result_rho1.2.rds"))
stopCluster(CL)

# rho 5 envir FALSE -------------------------------------------------------

library(parallel)
runs<- 1:1000
Nspp<- 50
Ncomm<- 50
power<- 5
perm<- 1000
envir<- FALSE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
})
ncor<- detectCores()
CL<- makeCluster(ncor, setup_timeout = 0.5)
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho5<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)
saveRDS(res_rho5, file = here::here("output", "Result_rho5.rds"))
stopCluster(CL)


# rho 0.0001 envir TRUE ------------------------------------------------------

runs<- 1:1000 # number of tests
Nspp<- 50 
Ncomm<- 50
power<- 0.0001 # grafen parameter
perm<- 1000 # number of permutations in taxa shuffle null model
envir<- TRUE # environment
u<- 5
binary<- TRUE
test<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
}) # simulating trees - equal to the number of tests

ncor<- detectCores() # setting number of cores
CL<- makeCluster(ncor, setup_timeout = 0.5) 
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

system.time(res_rho0.0001_envir<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)) # run analysis
stopCluster(CL)
saveRDS(res_rho0.0001_envir, file = here::here("output", "Result_rho0.0001_envir.rds"))



# rho 0.8 envir TRUE ------------------------------------------------------

runs<- 1:1000 # number of tests
Nspp<- 50 
Ncomm<- 50
power<- 0.8 # grafen parameter
perm<- 1000 # number of permutations in taxa shuffle null model
envir<- TRUE # environment
u<- 5
binary<- TRUE
test<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
}) # simulating trees - equal to the number of tests

ncor<- detectCores() # setting number of cores
CL<- makeCluster(ncor, setup_timeout = 0.5) 
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho.8_envir<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
) # run analysis
stopCluster(CL)
saveRDS(res_rho.8_envir, file = here::here("output", "Result_rho08_envir.rds"))

# rho 1.2 envir TRUE ------------------------------------------------------

library(parallel)
runs<- 1:1000
Nspp<- 50
Ncomm<- 50
power<- 1.2
perm<- 1000
envir<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
})
ncor<- detectCores()
CL<- makeCluster(ncor, setup_timeout = 0.5)
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho1.2_envir<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)
saveRDS(res_rho1.2_envir, file = here::here("output", "Result_rho1.2_envir.rds"))
stopCluster(CL)


# rho 5 envir TRUE --------------------------------------------------------

library(parallel)
runs<- 1:1000
Nspp<- 50
Ncomm<- 50
power<- 5
perm<- 1000
envir<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
})
ncor<- detectCores()
CL<- makeCluster(ncor, setup_timeout = 0.5)
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir", "perm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho5_envir<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = perm, power = power, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)
saveRDS(res_rho5_envir, file = here::here("output", "Result_rho5_envir.rds"))
stopCluster(CL)


