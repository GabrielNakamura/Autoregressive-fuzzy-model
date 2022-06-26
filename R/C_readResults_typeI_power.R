# load results from type I error and power analysis

# typeI error rho 0.0001 / no Env -------------------------------------------------

res_rho0001 <- readRDS(here::here("output", "Result_rho0001.rds"))
pval_0001 <- unlist(lapply(res_rho0001, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_0001 <= 0.05))/1000

# power rho 0.08 / no Env -------------------------------------------------

res_rho08<- readRDS(here::here("output", "Result_rho08.rds"))
pval_08<- unlist(lapply(res_rho08, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_08 <= 0.05))/1000

# power rho 1.2 / no Env -------------------------------------------------

res_rho1.2<- readRDS(here::here("output", "Result_rho1.2.rds"))
pval_1.2<- unlist(lapply(res_rho1.2, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_1.2 <= 0.05))/1000

# power rho 5 / no Env -------------------------------------------------

res_rho5<- readRDS(here::here("output", "Result_rho5.rds"))
pval_5<- unlist(lapply(res_rho5, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_5 <= 0.05))/1000

# typeI rho .0001 / Env -------------------------------------------------

res_rho0001_env<- readRDS(here::here("output", "Result_rho0.0001_envir.rds"))
pval_0001_env<- unlist(lapply(res_rho0001_env, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_0001_env <= 0.05))/1000
