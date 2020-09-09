# load results from type I error and power analysis

# typeI error rho 0.0001 / no Env -------------------------------------------------

res_rho0001<- readRDS(here::here("output", "Result_rho0001.rds"))
pval_0001<- unlist(lapply(res_rho0001, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_0001 <= 0.05))/1000


# power rho 0.08 / no Env -------------------------------------------------

res_rho08<- readRDS(here::here("output", "Result_rho08.rds"))
pval_08<- unlist(lapply(res_rho08, function(x){
  x$Model.Results[ncol(x$Model.Results)]
}))
length(which(pval_08 <= 0.05))/1000
