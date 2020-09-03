##### type I error and power analysis 

# data, functions and libraries -------------------------------------------
library(parallel)
source(here::here("functions", "autoregressive_GN_01-09.R"))


# power - 1 envir TRUE ----------------------------------------------------

runs<- 1:1000 # number of tests
Nspp<- 50 
Ncomm<- 50
power<- 1 # grafen parameter
perm<- 1000 # number of permutations in taxa shuffle null model
envir<- TRUE # environment

tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
}) # simulating trees - equal to the number of tests

ncor<- detectCores() # setting number of cores
CL<- makeCluster(ncor, setup_timeout = 0.5) 
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho1_envir<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = 10, power = 1, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
  ) # run analysis
stopCluster(CL)


# power 0.0001 envir TRUE -------------------------------------------------

library(parallel)
runs<- 1:1000
Nspp<- 50
Ncomm<- 50
power<- 0.0001
perm<- 1000
envir<- TRUE
tree_list<- lapply(runs, function(x){
  geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,t=30,extinct=FALSE)
})
ncor<- detectCores()
CL<- makeCluster(ncor, setup_timeout = 0.5)
clusterExport(CL,c("simul_metacomm_autoregress", "Nspp", "Ncomm", "power", "envir"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
clusterEvalQ(CL,library(QuantPsyc))

res_rho1_envir<- parLapply(cl = CL, X = tree_list, fun = function(i) 
  simul_metacomm_autoregress(tree = i, Nspp = Nspp, 
                             Ncomm = Ncomm, perm = 10, power = 1, u = 5, 
                             envir = envir, binary = TRUE, test = TRUE)
)
stopCluster(CL)

