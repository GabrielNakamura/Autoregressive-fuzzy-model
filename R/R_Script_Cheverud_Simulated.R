download.file(url = "https://github.com/GabrielNakamura/Autoregressive-fuzzy-model/archive/master.zip", destfile = "Autoregressive-fuzzy-model.zip")
unzip(zipfile = "Autoregressive-fuzzy-model.zip")

source("Simul_ARM_11fev22_LD.R")

#Simulation with E=P, diag=F:
sim.rho.0.0001.PhyOnly.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=FALSE,
                power=0.0001,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagF_test_n50r1000_FINAL.RData")
sim.rho.0.8.PhyOnly.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=FALSE,
                  power=0.8,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagF_test_n50r1000_FINAL.RData")
sim.rho.1.2.PhyOnly.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=FALSE,
                  power=1.2,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagF_test_n50r1000_FINAL.RData")
sim.rho.5.PhyOnly.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=FALSE,
                  power=5,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagF_test_n50r1000_FINAL.RData")

#Simulation with E=P, diag=T:
sim.rho.0.0001.PhyOnly.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=T,
      power=0.0001,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagT_test_n50r1000_FINAL.RData")
sim.rho.0.8.PhyOnly.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=T,
        power=0.8,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagT_test_n50r1000_FINAL.RData")
sim.rho.1.2.PhyOnly.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=T,
        power=1.2,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagT_test_n50r1000_FINAL.RData")
sim.rho.5.PhyOnly.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=FALSE,diag=T,
      power=5,binary=TRUE,runs=1000,test=TRUE,nperm=999,parallel=NULL)
save.image("sim.OnlyPhy_diagT_test_n50r1000_FINAL.RData")


#######################################################################
#Simulation with E≠P diag=F

sim.rho.0.0001.E.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
      diag=F,power=0.0001,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagF_n50r1000_FINAL.RData")
sim.rho.0.8.E.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
    diag=F,power=0.8,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagF_n50r1000_FINAL.RData")
sim.rho.1.2.E.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
            diag=F,power=1.2,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagF_n50r1000_FINAL.RData")
sim.rho.5.E.diagF<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
            diag=F,power=5,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagF_n50r1000_FINAL.RData")

#Simulation with E≠P diag=T

sim.rho.0.0001.E.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
        diag=T,power=0.0001,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagT_n50r1000_FINAL.RData")
sim.rho.0.8.E.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
            diag=T,power=0.8,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagT_n50r1000_FINAL.RData")
sim.rho.1.2.E.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
            diag=T,power=1.2,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagT_n50r1000_FINAL.RData")
sim.rho.5.E.diagT<-simul.comm(Ncomm=50,Nspp=50,u=5,envir=TRUE,
            diag=T,power=5,runs=1000,test=TRUE,nperm=999,parallel = NULL)
save.image("sim.WithEnv_test_diagT_n50r1000_FINAL.RData")


sim.rho.0.0001.PhyOnly.diagF$Model.Results
sim.rho.0.8.PhyOnly.diagF$Model.Results
sim.rho.1.2.PhyOnly.diagF$Model.Results
sim.rho.5.PhyOnly.diagF$Model.Results
sim.rho.0.0001.PhyOnly.diagT$Model.Results
sim.rho.0.8.PhyOnly.diagT$Model.Results
sim.rho.1.2.PhyOnly.diagT$Model.Results
sim.rho.5.PhyOnly.diagT$Model.Results
sim.rho.0.0001.E.diagF$Model.Results
sim.rho.0.8.E.diagF$Model.Results
sim.rho.1.2.E.diagF$Model.Results
sim.rho.5.E.diagF$Model.Results
sim.rho.0.0001.E.diagT$Model.Results
sim.rho.0.8.E.diagT$Model.Results
sim.rho.1.2.E.diagT$Model.Results
head(sim.rho.5.E.diagT$Model.Results)

Raw.res<-list(sim.rho.0.0001.PhyOnly.diagF$Model.Results,
              sim.rho.0.8.PhyOnly.diagF$Model.Results,
              sim.rho.1.2.PhyOnly.diagF$Model.Results,
              sim.rho.5.PhyOnly.diagF$Model.Results,
              sim.rho.0.0001.PhyOnly.diagT$Model.Results,
              sim.rho.0.8.PhyOnly.diagT$Model.Results,
              sim.rho.1.2.PhyOnly.diagT$Model.Results,
              sim.rho.5.PhyOnly.diagT$Model.Results,
              sim.rho.0.0001.E.diagF$Model.Results,
              sim.rho.0.8.E.diagF$Model.Results,
              sim.rho.1.2.E.diagF$Model.Results,
              sim.rho.5.E.diagF$Model.Results,
              sim.rho.0.0001.E.diagT$Model.Results,
              sim.rho.0.8.E.diagT$Model.Results,
              sim.rho.1.2.E.diagT$Model.Results,
              sim.rho.5.E.diagT$Model.Results)

Result.Table<-matrix(NA,16,5,
  dimnames=list(c("Phy.0001.ARM","Phy.08.ARM","Phy.12.ARM","Phy.5.ARM",
                  "Phy.0001.P","Phy.08.P","Phy.12.P","Phy.5.P",
                  "E.0001.ARM","E.08.ARM","E.12.ARM","E.5.ARM",
                  "E.0001.P","E.08.P","E.12.P","E.5.P"),
                c("R2.mod","Rho","Reject.Rate","R2.res","Rho.res")))
runs=1000
  for (i in 1:length(Raw.res)){
    Result.Table[i,c(1:2,4:5)]<-colMeans(Raw.res[[i]][,c(1,4,6,9)])
    Result.Table[i,3]<-sum(ifelse(Raw.res[[i]][,5]<= 0.05, yes=1, no=0))/runs
  }
Result.Table    

Beta.res<-cbind(sim.rho.0.0001.PhyOnly.diagF$Model.Results[,4],
              sim.rho.0.8.PhyOnly.diagF$Model.Results[,4],
              sim.rho.1.2.PhyOnly.diagF$Model.Results[,4],
              sim.rho.5.PhyOnly.diagF$Model.Results[,4],
              sim.rho.0.0001.PhyOnly.diagT$Model.Results[,4],
              sim.rho.0.8.PhyOnly.diagT$Model.Results[,4],
              sim.rho.1.2.PhyOnly.diagT$Model.Results[,4],
              sim.rho.5.PhyOnly.diagT$Model.Results[,4],
              sim.rho.0.0001.E.diagF$Model.Results[,4],
              sim.rho.0.8.E.diagF$Model.Results[,4],
              sim.rho.1.2.E.diagF$Model.Results[,4],
              sim.rho.5.E.diagF$Model.Results[,4],
              sim.rho.0.0001.E.diagT$Model.Results[,4],
              sim.rho.0.8.E.diagT$Model.Results[,4],
              sim.rho.1.2.E.diagT$Model.Results[,4],
              sim.rho.5.E.diagT$Model.Results[,4])    
write.table(Beta.res,"Beta_FINAL.txt",sep=" ")
