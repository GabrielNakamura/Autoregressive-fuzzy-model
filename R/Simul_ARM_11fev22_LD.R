simul.comm <- function(Ncomm,
                       Nspp,
                       envir=FALSE,
                       diag=FALSE,
                       u=5,
                       power=5,
                       binary=FALSE,
                       runs=30,
                       test=FALSE,
                       nperm=999,
                       parallel=8){
# Compute matrix P:
    matrix.p <- function(L, phylo){
        match<-picante::match.phylo.comm(phylo,L)
        matrix.w <- as.matrix(match$comm)
        similar.phy<-ape::vcv(phy=match$phy,corr=TRUE)
        matrix.phy <- 1/colSums(similar.phy)
        matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
          if(diag==FALSE){
            diag(matrix.q)<-0
          } else { matrix.q=matrix.q
            }
        matrix.P <- matrix.w %*% matrix.q
        return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
      }
    
# Double center matrix transformation: 
    matrix.double.center <- function(mat){
      mean_row_partial <- apply(mat, MARGIN = 2, 
                                function(x){
                                  x - rowMeans(mat)
                                }
      )
      mean_col_partial <- t(apply(mean_row_partial, MARGIN = 1, 
                                  function(x){
                                    x - colMeans(mat)
                                  }))
      double_center_matrix <- mean_col_partial + mean(as.matrix(mat))
      return(double_center_matrix)
    }

  #Res<-matrix(NA,nrow=runs,ncol=11,dimnames=list(1:runs,paste(c("R2.W","Coef.a.W","Coef.b.W","Beta.W","p.value","Value.PCPS.1","Value.PCPS.2","R2.Res","Coef.a.Res","Coef.b.Res","Beta.Res"))))
  Res<-matrix(NA,nrow=runs,ncol=9,dimnames=list(1:runs,paste(c("R2.W","Coef.a.W","Coef.b.W","Beta.W","p.value","R2.Res","Coef.a.Res","Coef.b.Res","Beta.Res"))))
  tree.list<-list()
  N<-matrix(NA,nrow=Nspp,ncol=runs,dimnames=list(1:Nspp,1:runs))
  K.N<-matrix(NA,nrow=runs,ncol=1,dimnames=list(1:runs,"K"))
  E<-matrix(NA,nrow=Ncomm,ncol=runs,dimnames=list(1:Ncomm,1:runs))
  L.list<-list()
  P.list<-list()
  pred.L.list<-list()
  resid.L.list<-list()
    for(k in 1:runs){
      phylo<-geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp,extinct=FALSE)
      tree.list[[k]]<-phylo
      a<-rnorm(Nspp,u,10)
      h<-runif(Nspp,0,30)
      L<-matrix(NA,Ncomm,Nspp,dimnames=list(paste("C", 1:Ncomm, sep=""),phylo$tip.label))
      rownames(N)<-phylo$tip.label
 	 
# Simulate niche N:    
      niche<-ape::rTraitCont(ape::compute.brlen(phylo,power=power),model="BM")
      N[,k]<-scales::rescale(niche,c(-1,101))
 	    K.N[k,]<-picante::Kcalc(N[,k],phylo)
# Simulate environment E:
 	    	 if(envir==TRUE){
	          E[,k]<-runif(Ncomm,0,100)
	       } else {E[,k]=N[,k]
	         }
# Simulate L:	    
		    for(i in 1:Ncomm){
			    for (j in 1:Nspp){
				    L[i,j]<-rpois(1,h[j]*exp((-((E[,k][i]-N[,k][j])^2))/(2*(a[j]^2))))
			    }
		    }
	    colnames(L)<-phylo$tip.label
	    rownames(L)<-rownames(as.matrix(L),FALSE,prefix="comm")
	      if(binary==TRUE){
	        L<-vegan::decostand(L,"pa")
	      } else{L=L}
	    L.list[[k]]<-L
	    L.cent<-matrix.double.center(L)

# Compute PCPS (silenced):
	    #pcps.L<-PCPS::pcps(L, cophenetic(phylo),checkdata = FALSE)
	    #Res[k,6]<-pcps.L$values[1,2]
	    #Res[k,7]<-pcps.L$values[2,2]

# Compute matrix P:	    
	    P<-matrix.p(L,phylo=phylo)$matrix.P
	    P.list[[k]]<-P
	    L.list[[k]]<-L
      P.cent<-matrix.double.center(P)

# run ARM:
      mod.L<-lm(as.numeric(L.cent)~as.numeric(P.cent))
	    Res[k,1]<-summary(mod.L)$r.squared
	    Res[k,2]<-summary(mod.L)$coefficients[1]
	    Res[k,3]<-summary(mod.L)$coefficients[2]
	    Res[k,4]<-QuantPsyc::lm.beta(mod.L)

# power/type I error test####     
	      if(test==TRUE){
	        seqpermutation.taxa <- SYNCSA::permut.vector(ncol(L), nset = nperm)
	        seqpermutation.taxa <- lapply(seq_len(nrow(seqpermutation.taxa)), function(i) seqpermutation.taxa[i,])
	        col.names<-colnames(L)
	        p.n.taxa <- function(samp, L, phylo){
	          L.null<-L
	          colnames(L.null)<-col.names[samp]
	          match.null<-picante::match.phylo.comm(phylo,L.null)
	          MP.null <- matrix.p(match.null$comm,match.null$phy)$matrix.P
	          return(MP.null)
	        }
	        P.null <- lapply(seqpermutation.taxa, p.n.taxa, L = L, phylo = phylo)
	        P.null.cent_list<- lapply(P.null, FUN=function(x){
	          P.null.cent<-matrix.double.center(x)
	        })

# start parallel:
	        newClusters <- FALSE
	        if (is.numeric(parallel)) {
	          parallel <- parallel::makeCluster(parallel, type = "PSOCK")
	          newClusters <- TRUE
	        }
	        if (!inherits(parallel, "cluster")) {
	          mod.L.null_list<- lapply(P.null.cent_list, function(x){
	            mod.L.null_list<- lm(as.numeric(L.cent) ~ as.numeric(x))
	          }
	          )
	        } else{
	          mod.L.null_list<- parallel::parLapply(cl = parallel,
	                                                X = P.null.cent_list, 
	                                                fun = function(i){
	                                                  lm(as.numeric(L.cent) ~ as.numeric(i))
	                                                }
	          )
	        }
	        if (newClusters){
	          parallel::stopCluster(parallel)
	        }
	        
	        mod.beta.null <- matrix(unlist(lapply(mod.L.null_list, function(x){
	                 QuantPsyc::lm.beta(x)})),
	                 nrow = nperm, 
	                 ncol= 1)
	        p.value<- (sum(ifelse(mod.beta.null > Res[k,4], 1, 0)) +1)/(nperm + 1)
	      } else {p.value=NA
	        }
	    Res[k,5]<-p.value
# Compute predicted and residual W matrices:
	    pred.L<-matrix(predict(mod.L),nrow(L),ncol(L),byrow=TRUE,dimnames=list(rownames(L),colnames(L)))
	    pred.L.list[[k]]<-pred.L
	    resid.L<-L.cent-pred.L
	    resid.L.list[[k]]<-resid.L
	    resid.L.cent<-matrix.double.center(resid.L)
	    mod.resid.L<-lm(as.numeric(resid.L.cent)~as.numeric(P.cent))
      Res[k,6]<-summary(mod.resid.L)$r.squared
      Res[k,7]<-summary(mod.resid.L)$coefficients[1]
      Res[k,8]<-summary(mod.resid.L)$coefficients[2]
      Res[k,9]<-QuantPsyc::lm.beta(mod.resid.L)
      print(k)
    }
  return(list(Trees=tree.list,Niche=N,K.niche=K.N,Environment=E,W.matrices=L.list,P.matrices=P.list,Predicted.W.matrices=pred.L.list,Residual.W.matrices=resid.L.list,Model.Results=Res))
  }